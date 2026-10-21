# Public-shell review fixes — 2026-10-21

Candidate and catalog detail buttons now use a dialog included in the shipped
local HTML. A Chromium regression loads that HTML and its real modules, opens
both detail paths, closes/reopens them, distinguishes input/candidate/reference/
redocking roles, and checks an external purchase link in a new tab. API records
and the supplier destination are synthetic; no supplier network request or
purchase is made by this test.

The table and new Excel exports explain the existing priority calculation:
the score from channels with values, the weight filled with neutral 50, the
saved calculation scope, and goal-dependent neutralization of the topical
term. These additions do not change numerical scores, ranking weights,
filters, pose QC, or saved studies. Raw values retain their full precision.

Stored product URLs appear under **구매처 홈페이지 열기** in candidate details.
New Excel exports include a clickable first purchase URL in the review table
and separate supplier links in each compound tab. Only stored HTTP(S) links
without embedded credentials are used. A link or catalog record does not
confirm current stock, price, purity, or shipping availability. Previous
downloaded workbooks and sealed archives are not rewritten.

Other regression fixes:

- Singleton ADMET objects and multi-model lists use the same readers; malformed
  items remain unassessed and original predictions are preserved.
- Frozen screening policies validate their own digest and supported semantics.
  Submission source hashes remain recorded; current execution hashes bind new
  cache keys and checkpoints. Changed execution/checkpoint conditions produce
  an explicit hold rather than silently mixing calculations.
- Legacy RCSB metadata retrieval sends at most 500 experimental IDs per call,
  merges validated batches and rejects incomplete/error responses. Its bounded
  regression uses 1,001 synthetic IDs, not a live 1,001-entry download.
- Large exact-name result sets can proceed after an explicitly selected CID is
  checked against the full returned identity set.
- Expired local-session errors explain refreshing the browser and checking
  saved studies before resubmitting a calculation.

Public tutorial screenshots use reviewed public-reference examples. Private
research results, databases and recipient delivery archives are excluded.

Verification uses isolated PostgreSQL, synthetic export fixtures and a real
Chromium DOM. Full-library screening, new native docking, GPU/OpenMM, OAuth
and biological validation are outside this repair run.

## Subsequent live link check

The reviewed purchase-link presentation was also applied to the existing
operating UI and exporter. One stored supplier destination was opened from
the live candidate detail in a new browser tab. The same destination was
opened using Microsoft Excel's **Open Hyperlink** command on an actual
exported candidate-table cell. This checked navigation to the supplier's
product page; it did not place an order or certify inventory.

The workbook check used a read-only copy of a completed operating snapshot
and the qualified exporter on a separate computer. This is distinct from
completing a new export job on the operating coordinator. Existing sealed
workbooks and archives are retained unchanged; purchase hyperlinks in a new
export do not retroactively add links to earlier files. Private inputs,
compound identifiers, outputs and screenshots are kept outside this repository.
