# Python 2 to Python 3 Port Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Migrate a Python 2 codebase to fully compatible, tested Python 3.

**Architecture:** Audit the codebase with automated tools, fix syntax and semantic incompatibilities in order of risk, update dependencies, then verify with a full test suite running under Python 3.

**Tech Stack:** Python 3.x, pytest, pyupgrade, 2to3 (stdlib), pip, optionally modernize/pylint

---

## Prerequisites

Before starting, ensure you have:
- Python 3.8+ installed (`python3 --version`)
- A virtual environment for Python 3 (`python3 -m venv .venv && source .venv/bin/activate`)
- The original test suite passing under Python 2 (baseline)

---

### Task 1: Audit — Identify All Incompatibilities

**Files:**
- Read: all `.py` files in the repo

**Step 1: Install audit tools**

```bash
pip install pyupgrade pylint
```

**Step 2: Run 2to3 in dry-run mode (no changes yet)**

```bash
2to3 --list-fixes .
2to3 -l  # list all available fixers
2to3 -v --no-diffs . 2>&1 | tee /tmp/py2to3-audit.txt
```

Expected: list of files and lines with incompatibilities. No files are modified.

**Step 3: Run pylint for Python 3 compatibility**

```bash
pylint --py3k **/*.py 2>&1 | tee /tmp/pylint-py3k.txt
```

**Step 4: Review and categorize findings**

Open `/tmp/py2to3-audit.txt` and `/tmp/pylint-py3k.txt`. Group issues into:
- Print statements
- Unicode/string handling
- Integer division
- Renamed/moved stdlib modules
- Dict `.keys()`, `.values()`, `.items()` misuse
- `xrange` usage
- Exception syntax
- `reduce`, `filter`, `map` behavior changes
- Dependency compatibility

**Step 5: Commit the audit output**

```bash
git add /tmp/py2to3-audit.txt /tmp/pylint-py3k.txt
git commit -m "chore: add Python 3 audit output"
```

---

### Task 2: Set Up Python 3 Test Environment

**Files:**
- Modify: `requirements.txt` or `setup.py` / `pyproject.toml`
- Create: `.python-version` (if using pyenv)

**Step 1: Write a failing smoke test**

```python
# tests/test_py3_compat.py
import sys

def test_python_version():
    assert sys.version_info.major == 3, "Must run under Python 3"
```

**Step 2: Run it to confirm it fails under Python 2**

```bash
python2 -m pytest tests/test_py3_compat.py -v
```

Expected: FAIL — "Must run under Python 3"

**Step 3: Run it to confirm it passes under Python 3**

```bash
python3 -m pytest tests/test_py3_compat.py -v
```

Expected: PASS

**Step 4: Check all dependencies for Python 3 support**

```bash
pip install pip-check-reqs
pip3 install -r requirements.txt 2>&1 | tee /tmp/dep-check.txt
```

Replace any Python 2-only packages (e.g. `mock` -> use `unittest.mock`, `futures` -> stdlib `concurrent.futures`).

**Step 5: Commit**

```bash
git add tests/test_py3_compat.py requirements.txt
git commit -m "chore: add Python 3 smoke test and check dependencies"
```

---

### Task 3: Fix Print Statements

**Files:**
- Modify: any `.py` file containing `print` without parentheses

**Step 1: Write a test that imports a module using print**

Identify the first module with print statements, e.g. `src/utils.py`.

```python
# tests/test_utils_import.py
def test_utils_importable():
    import src.utils  # will fail if print statement present
```

**Step 2: Run to confirm it fails**

```bash
python3 -m pytest tests/test_utils_import.py -v
```

Expected: SyntaxError

**Step 3: Apply the print fix with 2to3**

```bash
2to3 -w -f print src/
```

This converts `print x` to `print(x)` in all files under `src/`.

**Step 4: Run the test to confirm it passes**

```bash
python3 -m pytest tests/test_utils_import.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add src/
git commit -m "fix: convert print statements to print() calls"
```

---

### Task 4: Fix Unicode and String Handling

**Files:**
- Modify: any `.py` file using `unicode()`, `u""` prefix (in Python 2 style), or `str` vs `bytes` mismatches

**Step 1: Write a failing test for string behavior**

```python
# tests/test_strings.py
def test_str_is_unicode():
    s = "hello"
    assert isinstance(s, str)
    assert not isinstance(s, bytes)

def test_bytes_encoding():
    b = b"hello"
    assert isinstance(b, bytes)
    assert b.decode("utf-8") == "hello"
```

**Step 2: Run to confirm current behavior**

```bash
python3 -m pytest tests/test_strings.py -v
```

**Step 3: Remove `u""` prefix and `unicode()` calls**

```bash
2to3 -w -f unicode src/
pyupgrade --py3-plus src/**/*.py
```

Remove explicit `unicode()` calls — in Python 3, `str` is already unicode.

Replace `str(x)` with `x` where `x` is already a string, and `unicode(x)` with `str(x)`.

For any raw bytes handling, ensure explicit `.encode()` / `.decode()` is used.

**Step 4: Run tests**

```bash
python3 -m pytest tests/test_strings.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add src/
git commit -m "fix: update string handling for Python 3 unicode model"
```

---

### Task 5: Fix Integer Division

**Files:**
- Modify: any `.py` file using `/` division where integer truncation was assumed

**Step 1: Write a failing test**

```python
# tests/test_division.py
def test_integer_division():
    # In Python 3, / always returns float
    assert 5 / 2 == 2.5
    assert 5 // 2 == 2  # use // for integer division
```

**Step 2: Run to see current behavior**

```bash
python3 -m pytest tests/test_division.py -v
```

**Step 3: Apply fix with 2to3**

```bash
2to3 -w -f division src/
```

This adds `from __future__ import division` where needed, or you can manually replace `/` with `//` where integer division is intended.

**Step 4: Audit remaining `/` usages manually**

```bash
grep -rn " / " src/ --include="*.py"
```

For each hit, determine if float or integer division is intended and use `//` accordingly.

**Step 5: Run tests**

```bash
python3 -m pytest tests/test_division.py -v
```

Expected: PASS

**Step 6: Commit**

```bash
git add src/
git commit -m "fix: replace ambiguous / with // for integer division"
```

---

### Task 6: Fix Renamed Standard Library Modules

**Files:**
- Modify: any `.py` file importing renamed modules

**Step 1: Write failing import tests**

```python
# tests/test_stdlib_imports.py
def test_urllib_imports():
    from urllib.request import urlopen
    from urllib.parse import urlencode, urlparse

def test_configparser():
    import configparser

def test_queue():
    import queue

def test_http():
    import http.client
```

**Step 2: Run to see failures**

```bash
python3 -m pytest tests/test_stdlib_imports.py -v
```

**Step 3: Apply stdlib renames**

```bash
2to3 -w -f urllib -f urllib2 -f httplib -f ConfigParser -f Queue -f cPickle src/
```

Common renames:
| Python 2 | Python 3 |
|---|---|
| `urllib2` | `urllib.request`, `urllib.error` |
| `urllib.urlencode` | `urllib.parse.urlencode` |
| `httplib` | `http.client` |
| `ConfigParser` | `configparser` |
| `Queue` | `queue` |
| `cPickle` | `pickle` |
| `StringIO` | `io.StringIO` |
| `cStringIO` | `io.StringIO` |

**Step 4: Run tests**

```bash
python3 -m pytest tests/test_stdlib_imports.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add src/
git commit -m "fix: update stdlib imports to Python 3 module names"
```

---

### Task 7: Fix Dictionary Method Usage

**Files:**
- Modify: any `.py` file using `.keys()`, `.values()`, `.items()` where a list was expected

**Step 1: Write a failing test**

```python
# tests/test_dict_methods.py
def test_dict_keys_are_view():
    d = {"a": 1, "b": 2}
    keys = d.keys()
    assert "a" in keys  # views support 'in'
    keys_list = list(d.keys())  # explicit list conversion
    assert keys_list == ["a", "b"] or set(keys_list) == {"a", "b"}
```

**Step 2: Search for list-assumed dict method usage**

```bash
grep -rn "\.keys()\|\.values()\|\.items()" src/ --include="*.py"
```

Look for patterns like:
- `d.keys()[0]` — indexing a view (broken in Python 3)
- `len(d.keys())` — use `len(d)` instead
- Iterating is fine, but sorting requires `sorted(d.keys())`

**Step 3: Fix indexing and sorting**

Replace:
```python
d.keys()[0]          ->  list(d.keys())[0]
d.values()[0]        ->  list(d.values())[0]
sorted(d.keys())     ->  sorted(d)  # or sorted(d.keys())
```

**Step 4: Run tests**

```bash
python3 -m pytest tests/test_dict_methods.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add src/
git commit -m "fix: update dict method usage for Python 3 view objects"
```

---

### Task 8: Fix xrange, reduce, filter, map

**Files:**
- Modify: any `.py` file using `xrange`, `reduce`, `filter`, `map`

**Step 1: Write failing tests**

```python
# tests/test_builtins.py
def test_range():
    r = range(10)
    assert list(r) == list(range(10))

def test_filter_returns_iterator():
    result = list(filter(lambda x: x > 2, [1, 2, 3, 4]))
    assert result == [3, 4]

def test_map_returns_iterator():
    result = list(map(lambda x: x * 2, [1, 2, 3]))
    assert result == [2, 4, 6]

def test_reduce():
    from functools import reduce
    result = reduce(lambda a, b: a + b, [1, 2, 3])
    assert result == 6
```

**Step 2: Apply fixes**

```bash
2to3 -w -f xrange -f reduce -f filter -f map src/
```

- `xrange` -> `range`
- `reduce` -> `from functools import reduce`
- `filter(f, lst)` -> `list(filter(f, lst))` if result used as list
- `map(f, lst)` -> `list(map(f, lst))` if result used as list

**Step 3: Run tests**

```bash
python3 -m pytest tests/test_builtins.py -v
```

Expected: PASS

**Step 4: Commit**

```bash
git add src/
git commit -m "fix: replace xrange, reduce, update filter/map for Python 3"
```

---

### Task 9: Fix Exception Syntax

**Files:**
- Modify: any `.py` file using Python 2 exception syntax

**Step 1: Search for old exception syntax**

```bash
grep -rn "except.*," src/ --include="*.py"
```

Python 2: `except ValueError, e:`
Python 3: `except ValueError as e:`

**Step 2: Apply fix**

```bash
2to3 -w -f except src/
```

Also check for `raise Exception, "message"` syntax:
```bash
grep -rn "^raise " src/ --include="*.py"
```

Python 2: `raise ValueError, "message"`
Python 3: `raise ValueError("message")`

**Step 3: Write a test**

```python
# tests/test_exceptions.py
def test_exception_handling():
    try:
        raise ValueError("test")
    except ValueError as e:
        assert str(e) == "test"
```

**Step 4: Run tests**

```bash
python3 -m pytest tests/test_exceptions.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add src/
git commit -m "fix: update exception syntax to Python 3 style"
```

---

### Task 10: Run Full Test Suite and Fix Remaining Issues

**Files:**
- Modify: any remaining files with failures

**Step 1: Run the full test suite**

```bash
python3 -m pytest tests/ -v 2>&1 | tee /tmp/full-test-run.txt
```

**Step 2: Triage failures**

For each failing test, identify the root cause and fix it. Common remaining issues:
- `__future__` imports that are now no-ops (remove them)
- `has_key()` -> `in` operator
- `next()` vs `.next()` on iterators
- `zip()`, `enumerate()` returning iterators
- `input()` vs `raw_input()` (2to3 handles this)
- `long` type removed (use `int`)
- `file()` builtin removed (use `open()`)

Apply remaining 2to3 fixes:

```bash
2to3 -w -f has_key -f next -f zip -f input -f long -f file src/
```

**Step 3: Re-run full test suite**

```bash
python3 -m pytest tests/ -v
```

Expected: all tests PASS

**Step 4: Commit**

```bash
git add src/
git commit -m "fix: resolve remaining Python 3 incompatibilities"
```

---

### Task 11: Update Dependencies and Requirements

**Files:**
- Modify: `requirements.txt`, `setup.py`, or `pyproject.toml`

**Step 1: Write a test that imports all top-level dependencies**

```python
# tests/test_deps.py
def test_all_deps_importable():
    import importlib, pkg_resources
    for dist in pkg_resources.working_set:
        # just verify nothing crashes on import of the package
        pass
```

**Step 2: Pin Python 3 compatible versions**

Review each dependency in `requirements.txt`:
- Remove `mock` (use `unittest.mock`)
- Remove `futures` (use `concurrent.futures`)
- Remove `six` if it was only needed for Python 2/3 compatibility
- Update packages that have Python 3 versions (e.g. `Pillow` instead of `PIL`)

**Step 3: Install and test**

```bash
pip3 install -r requirements.txt
python3 -m pytest tests/ -v
```

Expected: PASS

**Step 4: Update setup.py / pyproject.toml classifiers**

```python
# In setup.py, update classifiers:
classifiers=[
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
],
python_requires=">=3.8",
```

**Step 5: Commit**

```bash
git add requirements.txt setup.py
git commit -m "chore: update dependencies and metadata for Python 3"
```

---

### Task 12: Final Verification

**Step 1: Run full test suite one last time**

```bash
python3 -m pytest tests/ -v --tb=short
```

Expected: all PASS, zero warnings about Python 2 syntax

**Step 2: Run pyupgrade for modern Python 3 style**

```bash
pyupgrade --py38-plus src/**/*.py
```

This removes legacy shims and upgrades to modern Python 3 idioms.

**Step 3: Run pylint**

```bash
pylint src/ --py3k
```

Expected: no Python 3 compatibility warnings

**Step 4: Final commit**

```bash
git add src/
git commit -m "chore: apply pyupgrade for modern Python 3 idioms"
```

**Step 5: Tag the release**

```bash
git tag -a v-py3 -m "Python 3 port complete"
```

---

## Common Pitfalls

| Issue | Symptom | Fix |
|---|---|---|
| Bytes vs str confusion | `TypeError: can't concat str to bytes` | Use `.encode()` / `.decode()` explicitly |
| Dict ordering assumed | Tests fail non-deterministically | Use `sorted()` or `collections.OrderedDict` |
| Integer division | Off-by-one or wrong results | Replace `/` with `//` where ints expected |
| `super()` calls | `TypeError: super() takes at least 1 argument` | Change `super(Class, self)` to `super()` |
| Metaclass syntax | `SyntaxError` | Change `__metaclass__ = X` to `class Foo(metaclass=X)` |

---

## References

- [Python 3 porting guide](https://docs.python.org/3/howto/pyporting.html)
- [2to3 fixer list](https://docs.python.org/3/library/2to3.html#fixers)
- [What's new in Python 3](https://docs.python.org/3/whatsnew/3.0.html)
