## R Standards

1. **Environment Management**
   - Always declare R dependencies in `DESCRIPTION` and lock versions with `renv.lock` for reproducibility.
   - Never hardcode secrets; use environment variables (`Sys.getenv()`) or secure key stores.

2. **Error Handling & Messaging**
   - Use `stop()` for errors, `warning()` for recoverable issues, and `message()` for informational logs.
   - Wrap external resources (DB connections, options, working directory changes) with `on.exit(..., add = TRUE)` for safe cleanup.

3. **Testing**
   - All new code must include unit tests with `testthat`.
   - Maintain >80% code coverage for critical transformations and SQL/data-mapping paths.

4. **Code Style**
   - Follow tidyverse style with `styler` and lint with `lintr`; keep functions small, explicit, and readable.
   - Prefer vectorized operations and type-stable outputs; use typed missing values (`NA_real_`, `NA_character_`) where appropriate.

5. **Documentation**
   - Document scripts/modules with clear headers and intent, inputs, outputs, and side effects.
   - Document exported/reusable functions with `roxygen2` tags (`@param`, `@return`, `@examples`).

6. **Data & SQL Safety**
   - Use parameterized queries via `DBI` (`dbGetQuery(..., params = ...)`); never build SQL with string interpolation for user/input values.
   - Preserve release/version guards for external data sources (e.g., Ensembl BioMart hosts) and avoid silent fallback behavior.

7. **CI/CD**
   - R workflows must run formatting/linting, tests, and parse checks on every PR.
   - For GitHub Actions in this repo, use shared ARC runners: `runs-on: ['atmos-aws-arc-runner-set']`.

8. **Dependency & Path Management**
   - Pin package versions with `renv`, review updates regularly, and validate compatibility after upgrades.
   - Avoid `setwd()` in reusable code; prefer project-relative paths and explicit path handling.
