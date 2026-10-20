# Contributing

Keep this repository limited to the SkinPipe platform: code, source-download metadata, empty configuration templates, documentation and synthetic tests.

Do not commit or attach research inputs, target lists from studies, candidate results, raw/processed chemical databases, predictions, private deployment settings, model weights, credentials, logs or collaborator delivery packages. Use temporary folders outside the checkout for runtime tests. Report a defect with a small synthetic reproducer and remove identifying data.

`PUBLIC_FILES.txt` is the reviewed source-file allowlist. Add a new software file deliberately to this list, then run:

```sh
python3 scripts/check-public-tree.py
python3 -m unittest discover -s tests -v
```

The guard checks the allowlist, empty study configuration and sensitive content patterns. It is an additional check, not a substitute for reviewing the diff. Build packages outside the checkout and run the package guard before distribution:

```sh
python3 scripts/check-public-tree.py --package /path/to/skinpipe.whl
```

Keep generated files and downloaded data outside Git. Do not merge the history of a research repository into this repository. Project distribution licensing remains undecided; retain all existing third-party notices.
