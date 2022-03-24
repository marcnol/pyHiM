# Good practices for commit

## Before
- Check if your code run without error
- Check if your code is clean (without debug `print()`, commented code line, ...)

## Content
- A commit should be **as atomic as possible**: 1 task/1 bug fix = 1 commit
- A commit must **always contain a message** explaining the reason and/or the content of the commit

## Message

To improve readability, commit message should follow this template:

```
<context>(<scope>) - #<issue_number> Summary title

- List of task or bugfix done in this commit
- Each item must start with a markdown bullet character like `-`
```

- **<context>**:
    - `fix` for bug fix
    - `feature` for development about a new feature
    - `doc` for adding documentation
    - `refactor` for improvement of readability, reusability or structure

- **<scope>** is the feature ou field involved by this commit (optionnal)

- `<issue_number>` refer to GitHub issues (optionnal)

**Example:**
```text
doc(contributor) - #167 Write file for commit good practices

- Before check-list
- Content reminder
- Message template
- This example
```