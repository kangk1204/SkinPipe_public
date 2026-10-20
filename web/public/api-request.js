// Allow the server's bounded CLI work to finish before the browser deadline.
export function requestErrorText(error) {
  if(error?.name==='TimeoutError')return '서버 응답 시간이 초과되었습니다. 서버 연결을 확인하고, 실행을 요청했다면 작업·계산 자원에서 접수 상태를 먼저 확인하세요.';
  if(error?.name==='AbortError')return '요청이 중단되었습니다. 실행을 요청했다면 접수 상태를 먼저 확인하세요.';
  if(error?.name==='TypeError'&&/fetch|network|load failed/i.test(error.message||''))return '서버에 연결할 수 없습니다. 서버 연결을 확인해 주세요. 저장된 후보 결과는 계속 확인할 수 있습니다.';
  return error?.message||String(error);
}
export function requestTimeout(path) {
  // Account-token refresh and the bounded provider catalog request precede this reply.
  if (path === '/v2/literature/models') return 90000;
  if (/^\/v2\/exports\/[^/]+\/download$/.test(path)) return 1800000;
  if (path.startsWith('/v2/experiment-batches') || path === '/v2/plans' || path === '/v2/feedback-snapshots' || path === '/v2/catalog-snapshots' || path.startsWith('/v2/campaigns')) return 120000;
  if (path.startsWith('/v2/catalog-')) return 40000;
  if (path === '/v2/correspondence') return 40000;
  // A fresh server verified the complete 60-analysis index in 70 seconds.
  // Result reads keep all integrity checks and a separate, bounded IO budget.
  if (path === '/v1/screens' || path.startsWith('/v1/screens/')) return 120000;
  // First v2 identity lookup verifies the pinned local ChEMBL database bytes.
  if (path === '/v2/targets/resolve' || path === '/v2/research-inputs/resolve') return 600000;
  // Synchronous receptor downloads and validation precede the acceptance reply.
  if (path === '/v2/discovery-runs') return 600000;
  if (path === '/v2/research-inputs/uploads' || path === '/v2/research-assistance' || path === '/v2/literature/reviews') return 120000;
  if (path === '/v1/resolve-target' || path === '/v1/assess') return 100000;
  if (path === '/v1/compare-hops') return 40000;
  if (path === '/v1/hop-jobs') return 100000; // Original and parent integrity checks precede enqueue.
  return 20000;
}
export async function requestJson(fetcher, url, path, request) {
  const method = (request.method || 'GET').toUpperCase();
  const repeatable = method === 'GET' || (method === 'POST' &&
    ['/v1/resolve-target','/v1/compare-hops','/v2/candidate-queries'].includes(path));
  for (let attempt = 0; ; attempt++) {
    try {
      const response = await fetcher(url, request);
      if (!response.ok) {
        let message;
        try { message = (await response.json()).error; } catch {}
        const error=Error(message || `서버 응답 오류 (${response.status})`);
        error.status=response.status;
        throw error;
      }
      // A connection can fail after headers arrive, while the body is read.
      return await response.json();
    } catch (error) {
      if (attempt || !repeatable || error.name !== 'TypeError' || request.signal?.aborted) throw error;
    }
  }
}
