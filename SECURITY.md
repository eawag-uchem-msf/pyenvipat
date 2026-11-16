# Security Policy

## Supported versions

Only the latest released version of pyEnvipat receives security fixes.

## Reporting a vulnerability

Please do **not** report security issues in public issues or MRs.

Report privately by email to the maintainer, johannesboog@yahoo.de, or by
opening a GitLab issue marked **confidential**. Include the affected version,
a description, and a minimal reproduction (for example, an input formula that
triggers a crash). pyEnvipat contains a compiled C extension that parses
user-supplied formula strings, so memory-safety reports (out-of-bounds access,
crashes on malformed input) are in scope.

You can expect an acknowledgement within a few working days. We will
coordinate a fix and a release with you and credit you in the release notes
unless you prefer otherwise.

Vulnerabilities in the original R package or its upstream C sources should
also be reported to the upstream project, <https://github.com/blosloos/enviPat>.
