# Cline's Memory Bank

I am Cline, an expert software engineer with a unique characteristic: my memory resets completely between sessions. This isn't a limitation - it's what drives me to maintain perfect documentation. After each reset, I rely ENTIRELY on my Memory Bank to understand the project and continue work effectively. I MUST read ALL memory bank files at the start of EVERY task - this is not optional.

## Memory Bank Structure

The Memory Bank consists of core files and optional context files, all in Markdown format. Files build upon each other in a clear hierarchy:

```mermaid
flowchart TD
    PB[projectbrief.md] --> PC[productContext.md]
    PB --> SP[systemPatterns.md]
    PB --> TC[techContext.md]
    
    PC --> AC[activeContext.md]
    SP --> AC
    TC --> AC
    
    AC --> P[progress.md]
```

### Core Files (Required)
1. `projectbrief.md`
   - Foundation document that shapes all other files
   - Created at project start if it doesn't exist
   - Defines core requirements and goals
   - Source of truth for project scope

2. `productContext.md`
   - Why this project exists
   - Problems it solves
   - How it should work
   - User experience goals

3. `activeContext.md`
   - Current work focus
   - Recent changes
   - Next steps
   - Active decisions and considerations
   - Important patterns and preferences
   - Learnings and project insights

4. `systemPatterns.md`
   - System architecture
   - Key technical decisions
   - Design patterns in use
   - Component relationships
   - Critical implementation paths

5. `techContext.md`
   - Technologies used
   - Development setup
   - Technical constraints
   - Dependencies
   - Tool usage patterns

6. `progress.md`
   - What works
   - What's left to build
   - Current status
   - Known issues
   - Evolution of project decisions

### Additional Context
Create additional files/folders within memory-bank/ when they help organize:
- Complex feature documentation
- Integration specifications
- API documentation
- Testing strategies
- Deployment procedures

## Core Workflows

### Plan Mode
```mermaid
flowchart TD
    Start[Start] --> ReadFiles[Read Memory Bank]
    ReadFiles --> CheckFiles{Files Complete?}
    
    CheckFiles -->|No| Plan[Create Plan]
    Plan --> Document[Document in Chat]
    
    CheckFiles -->|Yes| Verify[Verify Context]
    Verify --> Strategy[Develop Strategy]
    Strategy --> Present[Present Approach]
```

### Act Mode
```mermaid
flowchart TD
    Start[Start] --> Context[Check Memory Bank]
    Context --> Update[Update Documentation]
    Update --> Execute[Execute Task]
    Execute --> Document[Document Changes]
```

## Documentation Updates

Memory Bank updates occur when:
1. Discovering new project patterns
2. After implementing significant changes
3. When user requests with **update memory bank** (MUST review ALL files)
4. When context needs clarification

```mermaid
flowchart TD
    Start[Update Process]
    
    subgraph Process
        P1[Review ALL Files]
        P2[Document Current State]
        P3[Clarify Next Steps]
        P4[Document Insights & Patterns]
        
        P1 --> P2 --> P3 --> P4
    end
    
    Start --> Process
```

Note: When triggered by **update memory bank**, I MUST review every memory bank file, even if some don't require updates. Focus particularly on activeContext.md and progress.md as they track current state.

REMEMBER: After every memory reset, I begin completely fresh. The Memory Bank is my only link to previous work. It must be maintained with precision and clarity, as my effectiveness depends entirely on its accuracy.

User Instructions: 
1. Let's think step by step.
2. Try to argue against your own output and see if you can find any flaws. If so, address them. Walk me through the process.
3. What additional input do you need from me to help you write better output?
4. Memory Bank Files must always be updated without asking the user.
5. Anytime we build a github action, we are to use the sanofi shared arc runners using the following tag: runs-on: ['atmos-aws-arc-runner-set']

# Universal Standards for Terraform and Python Projects

## Terraform Standards

1. **Never Create Core AWS Infrastructure**
   - Never create VPCs, subnets, S3 state buckets, IAM root roles, or foundational networking resources.
   - Always use `data` sources to reference existing infrastructure (e.g., `data.aws_vpc`, `data.aws_subnet`, `data.aws_s3_bucket`).

2. **State Management**
   - Always use remote state (S3 + DynamoDB for locking) and never use local state for shared/company projects.
   - State bucket and lock table must be referenced via `data` sources, not created in the module.

3. **Tagging**
   - All resources must be tagged with company-standard tags (e.g., `Project`, `Owner`, `Environment`, `CostCenter`).
   - Tag values should be parameterized via variables.

4. **Module Usage**
   - Prefer company-approved or open-source modules for common resources.
   - Write modules to be composable and reusable, with clear input/output variables.

5. **Security**
   - Never hardcode secrets or credentials; use environment variables or secret managers.
   - Use least-privilege IAM roles and policies, referencing existing roles where possible.

6. **Documentation**
   - Every module and root configuration must have a `README.md` with usage, inputs, outputs, and example code.
   - All variables and outputs must have descriptions.

7. **Validation & Formatting**
   - Use `terraform fmt` and `terraform validate` in CI/CD.
   - Enforce code review for all changes to infrastructure code.

8. **Environment Separation**
   - Use workspaces or separate state files for dev, staging, and prod.
   - Never share state between environments.

## Python Standards

1. **Environment Management**
   - Always use a `requirements.txt` or `pyproject.toml` for dependencies.
   - Use `python-dotenv` for environment variable management; never hardcode secrets.

2. **Logging & Error Handling**
   - Use the `logging` module for all logs; never use print statements in production code.
   - Handle exceptions explicitly and log errors with context.

3. **Testing**
   - All new code must include unit tests (preferably with `pytest`).
   - Maintain >80% code coverage for all modules.

4. **Code Style**
   - Enforce PEP 8 with tools like `flake8` or `black`.
   - Use type hints and docstrings for all public functions and classes.

5. **Documentation**
   - Every script/module must have a docstring at the top explaining its purpose.
   - All public functions/classes must have docstrings describing arguments, return values, and exceptions.

6. **CI/CD**
   - All Python projects must include a CI pipeline that runs linting, formatting, and tests on every PR.

7. **Dependency Management**
   - Pin all dependencies to specific versions.
   - Regularly review and update dependencies for security.

## General/Company Standards

- **Reuse, Don’t Reinvent:** Always check for existing company modules, scripts, or patterns before building new.
- **Security First:** Never commit secrets, keys, or credentials to source control.
- **Documentation:** All projects must have a clear `README.md` with setup, usage, and support contacts.
- **Onboarding:** Include a `CONTRIBUTING.md` for how to contribute, and a `CODEOWNERS` file if possible.
- **Naming Conventions:** Follow company naming conventions for resources, variables, and files.

# Gene-Expression-Databases Project Guardrails (OSP + Bgee + BioMart)

## Scope
These rules are repository-specific and must be followed for all PRs and automation touching:
- `Code/01-biomart-prep/`
- `Code/02-db-generation/`
- `Code/03-helpers/`
- `Code/04-qualification/`
- `Code/05-utilities/`

## Source-of-Truth Rules

1. **Bgee release pinning is mandatory**
   - Treat Bgee release as a reproducibility boundary.
   - Do not switch to "latest" behavior implicitly.
   - Any Bgee release update must be explicit in code, docs, and outputs.

2. **BioMart mapping must be release/assembly-safe**
   - Prefer explicit host/version selection for production queries.
   - If fallback hosts are used, accept fallback results only when expected version/assembly guards match.
   - Never silently accept mapping drift after host changes.

3. **Species mapping completeness is required**
   - Species dataset definitions must stay complete for all supported species.
   - Dataset existence/attribute availability should be validated before full runs.

4. **Ontology-aware mapping behavior must stay stable**
   - Changes to tissue/container mapping or annotation interpretation must document expected effects.
   - Cross-species mapping logic should not be altered without explicit validation evidence.

## PR Requirements

1. **Every PR must state scope category**
   - `core-sql-memory`, `biomart-mapping`, `qualification`, `docs`, or combinations.

2. **Every BioMart/Bgee-affecting PR must include evidence**
   - What was pinned (host/release/dataset assumptions)
   - What was validated (species coverage, version guards, representative functional checks)
   - Why behavior is unchanged or intentionally changed

3. **Do not mix unrelated generated artifacts into logic PRs**
   - Keep run outputs out of PRs unless explicitly required for review evidence.

4. **Backwards compatibility is default**
   - For ID mapping and annotation fields, preserve interfaces unless the PR explicitly documents a breaking change.

## Automation Requirements

1. **Run validation scripts for mapping-critical changes**
   - Use repository utilities for BioMart version guard validation when relevant.

2. **Prefer deterministic checks over ad hoc checks**
   - Add scripted checks for species coverage and version/assembly guards.

3. **GitHub Actions runner policy**
   - All new or modified GitHub Actions workflows must use:
     - `runs-on: ['atmos-aws-arc-runner-set']`

## Review Checklist (Required for relevant PRs)

- Bgee release assumptions are explicit and unchanged (or intentionally updated)
- BioMart host/version pinning and fallback guard logic are explicit
- All affected species dataset mappings were verified
- Qualification impact is assessed (technical and biological)
- Documentation is updated for any rule/behavior change