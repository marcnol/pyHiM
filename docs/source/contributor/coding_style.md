# Python coding style

- If you have any doubts, please refer to the [PEP8](https://www.python.org/dev/peps/pep-0008/).
- To apply an auto format for PEP8, use `black` linter.

## Some examples of [PEP8](https://www.python.org/dev/peps/pep-0008/)

### Names

**TO AVOID**, never use these characters as single character variable names : 

- 'l' (lowercase letter 'L' [el])
- 'O' (uppercase letter 'o' [oh])
- 'I' (uppercase letter 'i' [eye])

| Type                        | Style                                                |
| --------------------------- | ---------------------------------------------------- |
| module                      | `lowercase_with_underscores` && short !              |
| class name                  | `CapWords`                                           |
| Variable Names              | `lowercase_with_underscores` && the more explicit !! |
| Function and Variable Names | `lowercase_with_underscores` && start with verb      |
| Constant                    | `UPPERCASE_WITH_UNDERSCORES`                         |

### import order

1. Standard library imports.
2. Related third party imports.
3. Local application/library specific imports.

```python
import mypkg.sibling
from mypkg import sibling
from mypkg.sibling import example
```

## The Zen of Python ([PEP20](https://www.python.org/dev/peps/pep-0020/#id2))

*If you're looking for the Python spirit...*

```
Beautiful is better than ugly.
Explicit is better than implicit.
Simple is better than complex.
Complex is better than complicated.
Flat is better than nested.
Sparse is better than dense.
Readability counts.
Special cases aren't special enough to break the rules.
Although practicality beats purity.
Errors should never pass silently.
Unless explicitly silenced.
In the face of ambiguity, refuse the temptation to guess.
There should be one-- and preferably only one --obvious way to do it.
Although that way may not be obvious at first unless you're Dutch.
Now is better than never.
Although never is often better than *right* now.
If the implementation is hard to explain, it's a bad idea.
If the implementation is easy to explain, it may be a good idea.
Namespaces are one honking great idea -- let's do more of those!
```