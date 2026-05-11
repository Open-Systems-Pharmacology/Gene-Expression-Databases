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
