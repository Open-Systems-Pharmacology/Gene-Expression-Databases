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