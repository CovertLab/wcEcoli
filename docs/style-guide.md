# Style Guide

_These guidelines are condensed from the PEP8 and Google style guides and our team decisions. -- Jerry_

* See [PEP8](https://www.python.org/dev/peps/pep-0008/)
* See [Google's Python style guide](https://google.github.io/styleguide/pyguide.html)
* See our [March 12, 2018 meeting slides](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_55) with notes on the tradeoffs between alternatives.


**What's a style guide?** Recommendations on coding alternatives like naming style, code layout,
docstring format, and common code patterns.

**Purpose:**
* Avoid some classes of bugs.
* Increase code readability, familiarity, and portability.
* Reduce work when merging code changes.
* Make it easier to collaborate.

But don't overdo consistency. There are cases for exceptions.

This doc (and the 2018 meeting slides) should summarize the team's decisions on
each PEP8 guideline to:
* Expected (point it out in code reviews), or
* Soft target (don't sweat it in code reviews), or
* Not adopted.

and whether to:
* Adopt it for new code.
* Plan to change existing code gradually or rapidly.


## wcEcoli Style Guidelines from PEP8 ... with adjustments

### Misc

* Use **UTF-8 text encoding** (otherwise the source file needs an encoding declaration).

* Module content order:

  ```python
  """Module docstring."""
  __future__ imports
  __all__ = ['a', 'b', 'c']  # and other "dunder" settings like __version__ and __author__
  imports
  code
  ```

* Aim for **docstrings** on all public modules, functions, classes, and methods.
  * `"""Docstrings should use triple quotes."""` whether it's one line or many; `"""` or `'''`.
    (The `"""` that ends a multiline docstring is supposed to go on a line by itself for some
    reason.)
  * A function docstring should start with an imperative sentence like "Plot all the things.".
  * Add docstrings to existing classes and functions while working on the code.
    A short introduction can be a big help over nothing.
  * Comments and docstrings that contradict the code are worse than no comments.
  * Write complete sentences.
  * Tip: If a docstring contains LaTeX code, make it a raw docstring:

    ```python
    r"""Fits the values of alpha and r such that the computed RNA synthesis
    probabilities converge to the measured RNA synthesis probabilities.
    
    v_{synth, j} = \alpha_j + \sum_{i} P_{T,i}*r_{ij}
    """
    ```

* Comments are usually complete sentences. The first word should be capitalized unless it's
  an identifier that begins with a lower case letter.


### Indentation

* Stick with TABs in this project.
* Set your editor's TAB stops to 4 spaces.
* Python assumes TAB stops are at 8 spaces unless told otherwise.
* Python rejects mixing tab and space indentation in a file.
* PEP8 recommends 4 spaces per indentation level for _shared code_.  
  (In years of early experience, TABs caused endless problems for [shared code](http://bugs.python.org/issue7012#msg93280). -- Tim Peters)  
  (Thou shalt indent with four spaces. No more, no less. Four shall be the number of spaces thou shalt indent, and the number of thy indenting shall be four. Eight shalt thou not indent, nor either indent thou two, excepting that thou then proceed to four. Tabs are right out. -- Georg Brandl)


### Imports

* **Imports at the top of the file.** Definitely don't put `import` in repeatedly-executed code.

  Occasionally there are reasons to break this rule, like `import pdb` or a
  slow-loading import in rarely-used code.

* Import separate modules on separate lines.

* **Import order:** standard library imports, blank line, third party imports, blank line, local imports.  
   Alphabetizing each group by module name is a bit nicer and reduces code merge work.

* **Avoid wildcard imports** (`from <module> import *`) since those make it unclear what names
  are present in the namespace, confusing readers and some tools. It's reasonable
  when republishing an internal interface as part of a public API (for example, overwriting a
  pure Python implementation of an interface with the definitions from an optional accelerator
  module and exactly which definitions will be overwritten isn’t known in advance).

* There are recommendations to import only modules rather than names from modules to
  avoid certain kinds of problems. The wcEcoli project doesn't follow this recommendation
  but it might be sensible to move in that direction.
  It's conventional and convenient to import names from `typing`
  (`from typing import Any, cast`), `collections.abc`, and some other modules,
  and a case could be made to treat classes as scopes like modules in this sense.

* Mostly use **absolute imports**. Occasionally it's worth using an explicit relative import,
  esp. when an absolute import is unnecessarily verbose:

      from . import sibling
      from path.to.mypkg import sibling
      from .sibling.submodule import example


### Lines

* **Line length:** soft target at 79 columns; harder target at 99 columns; no hard limit. Ditto for comments and docstrings (unlike PEP8).
  * A standard line length aids editor window layout and diff displays, but bio simulations might have many long names. It's annoying to stay within a hard limit but useful to have a shared target.
  * (PEP8 recommends 79 columns for code and 72 for comments and docstrings.)

* Prefer Python's **implied line continuation** inside parentheses, brackets and braces over a
  backslash for line continuation. Backslashes get broken by invisible white space at the end
  of a line.
* Avoid trailing whitespace on a line -- it breaks backslash line continutations.

* Prefer this indentation for function parameters:

      def long_function_name(
              var_one, var_two, var_three,
              var_four):
          print(var_one)

  or maybe put the closing parenthesis on a separate line:

      def long_function_name(
              var_one, var_two, var_three,
              var_four
              ):
          print(var_one)

  or the PEP8 alternative:

      foo = long_function_name(var_one, var_two,
                               var_three, var_four)

  not:

      def long_function_name(var_one, var_two,
              var_three, var_four):
          print(var_one)

  PyCharm is configurable but it implements the first indentation style by default, and using
  its Refactor command to rename `long_function_name()` will re-indent the continuation lines.

* One blank line between method definitions. (PEP8 calls for two blank lines between top-level
  definitions but we aren't adhering to that.) An additional blank line is OK to separate sections.

* Prefer a line break before a binary operator rather than after.
  This math tradition is clearer:
  ```python
  income = (gross_wages
          + taxable_interest
          + (dividends - qualified_dividends)
          - ira_deduction
          - student_loan_interest)
  ```
* Do likewise with whitespace between lines of a multiline string.
  ```python
  self.define_option(parser, 'timestamp', str, fp.timestamp(),
      help='Timestamp for this workflow. It gets combined with the'
           ' Workflow ID to form the workflow name. Set this if you want'
           ' to upload new steps for an existing workflow. Default ='
           ' the current local date-time.')
  ```


### Spaces

* Put at least one space before an inline comment, preferably two, and one after: `␣␣#␣`. (PEP8 says at least two spaces before the `#`.)

* Spaces in expressions (see [PEP8](https://www.python.org/dev/peps/pep-0008/) for more info):

  ```python
  # No spaces immediately within `()`, `[]`, or `{}`.
  spam(ham[1], {eggs: 2, salsa: 10})
  
  # A space between `,` `;` or `:` and any following item.
  demo = (0,) + (2, 3)
  if x == 4:
      print x, y; x, y = y, x
  
  # No space in a simple slice expression, but parentheses to clarify
  # complicated slice precedence, or construct a slice object, or put
  # subexpressions into variables.
  ham[1:9], ham[1:9:3], ham[:9:3], ham[1::3], ham[1:9:]
  ham[(lower+offset):(upper+offset)]
  
  # Put no space in function application or object indexing.
  spam(1) < spam(2)
  dct['key'] += lst[index]

  # Don't line up the `=` on multiple lines of assignment statements.
  # Maintaining it costs work and merge conflicts.
  x = 1
  long_variable = (3, 10)
  
  # Spaces around kwarg `=` are OK, unlike in PEP8, which recommends them only
  # with a type annotation.
  c = magic(real=1.0, imag=10.5)
  c = magic(real = 1.0, imag = 10.5)
  def munge(input: AnyStr, sep: AnyStr = None, limit=1000): ...
  
  # Use spaces or parentheses to help convey precedence.
  # Put zero or one space on each side of a binary operator (except indentation).
  hypot2 = x*x + y*y

  x: int = 333
  ```

* Avoid compound statements on one line.

  ```python
  if foo == 'blah': do_something()
  ```


### Naming style

   ```python
   ClassName
   ExceptionName  # usually ends with "Error"
   
   GLOBAL_CONSTANT_NAME
   
   function_name, method_name
   decorator_name
   
   local_var_name, global_var_name, instance_var_name, function_parameter_name
   camelCase  # it's OK to match the existing style but it's good to change it to snake case

   _internal_name
   __mangled_class_attribute_name  # https://peps.python.org/pep-0008/#designing-for-inheritance
   
   module_name
   package  # underscores are discouraged
   ```

* Use a trailing `_` to avoid a name conflict with a Python keyword or global function e.g.
  `id_`, `yield_`, `complex_`, or `max_`.
* Public names follow the convention for how it's used rather than how it's implemented, e.g.
  a class used as a decorator or like a callable function.
* Always use `self` for the first argument of an instance method and `cls` for the first
  argument of a class method.
* Don't use `l`, `O`, or `I` for a single letter variable name because they're hard to
  distinguish in some fonts.
* Don't invent `__dunder__` names.
* It's recommended to not make exceptions for scientific naming like `Kcat` and math
  naming like matrix `M`, but wcEcoli doesn't stick to this recommendation.
* Avoid using properties for expensive operations. The notation `obj.property` suggests
  it's low cost. For expensive operations, use methods, say, `obj.get_something()`
  or `obj.compute_something()`.


## From Google's Python Style Guide

See https://google.github.io/styleguide/pyguide.html

* Don't put code at the top level that you don't want to run when the file is imported, unit
  tested, or pydoc'd.
* Avoid mutable global state (including global variables), metaclasses,
  bytecode access, dynamic inheritance, object reparenting, reflection, modifying
  system internals, and `__del__` methods implementing customized cleanup (except
  standard library features like `enum` that use them internally),
  (at least without very compelling reasons).
* Don't use mutable objects as default values of a function or method.
* Import only packages and modules, not names from modules.
* Run the pylint and/or PyCharm inspections.
* Use full pathnames to import a module. Avoid relative imports to help prevent importing a package twice.
* Prefer the operator module e.g. `operator.mul` over `lambda x, y: x * y`. There's also `operator.itemgetter(*items)`, `operator.attrgetter(*attrs)`, and `operator.methodcaller(name[, args...])`.
* Use Python falsey tests for empty sequences and 0, e.g. `if sequence:` rather than `if len(sequence):` or `if len(sequence) > 0:`, but not for testing if a value is (not) `None`.
  * But don't write `if value:` to test for a non-empty string. That can be confusing.
* Don't use parentheses in return statements or conditional statements except for implied line continuation.
  * PyCharm recommends `return x, y` over `return (x, y)`.

* Write a docstring as a summary sentence, a blank line, then the rest.
* A function must have a docstring, unless it's: not externally visible, short, and obvious.
  Explain what you need to know to call it.
* Classes should have docstrings.
* `TODO: Investigate optimizations`. Searchable "TODO" comment. Mention an Issue number if any.
  (The older style was `TODO(username):`, naming who to ask for more information.)
* Use `assert` in unit tests and to check internal correctness. For a mistake at the API level,
  raise an exception. Raise built-in exception classes like `ValueError` when it makes sense.


## Pythonic Style

```python
if x is None: ...  # better than `if x == None:`
if x is not None: ...
def f(x): return 2 * x  # allows type hints & better stack traces than `f = lambda x: 2 * x
isinstance(obj, int)  # better than `type(obj) is int`
' '.join(a_sequence)  # faster than a loop of incremental string concatenation
if sequence: ...  # preferred over `if len(sequence):`
def f(x=None):
  x = x or []   # solid, unlike a mutable default argument `def f(x=[]):`
  return x, x * 2, x * 3
```

* Use `if x is None:` or `if x is not None:` rather than `==` or `!=`. Ditto for
  enums and other singletons. It's faster, esp. if it avoids calling a custom
  `__eq__()` method, and it might avoid exceptions or incorrect results in `__eq__(None)`.
* Avoid manually testing data types. Prefer object-oriented dispatch or `functools`
  features like `singledispatch`, or at least `isinstance()` since it handles subtypes
  and protocols.
* All or none of a function's `return` statements should return a value, usually of a
  consistent type. You can use a `dict` or a class instance to handle a union of cases.
* Prefer `str` methods over the `string` module. They're faster.
* Use `' '.join()` rather than looping over `a_string += stuff` to combine strings since
  `join()` takes linear time whereas repeated string concatenation can take O(_n_^2) time
  (unless a particular CPython optimization kicks in).

* When implementing ordering operations with rich comparisons, it's best to implement all
  six operations or use `@functools.total_ordering()` to fill them out.
* Use `def f(x): return 2 * x` instead of `f = lambda x: 2 * x` to get more helpful
  stack traces and allow type hints.
* Prefer f-string or `str.format(...)` over printf-style `%`-formatting because
  [`%`-formatting has quirks](https://docs.python.org/3/library/stdtypes.html?highlight=sprintf#printf-style-string-formatting)
  that lead to common errors such as failing to display tuples and dictionaries
  correctly.

### Exceptions
* Use the bare `except:` clause only when just printing/logging the traceback and
  re-raising the exception, or in a debugging tool.
* Limit a `try` clause to a narrow range of code so it only doesn't bury totally
  unexpected exceptions.
* Use a `with` statement or `try`/`finally` to ensure cleanup gets done, like closing a file.
* Any kind of failure should raise an explicit exception.
* Derive exceptions from `Exception` rather than `BaseException` unless catching this
  exception is almost always the wrong thing to do, like `KeyboardInterrupt` and `SystemExit`.
* When raising exceptions, aim to answer "What went wrong?" rather than just
  indicating "A problem occurred."


# Type hints

Python type hints can help catch bugs, improve code documentation,
and aid dev tools such as PyCharm's code completion and refactoring.
They can also provide input to tools like RPC parameter marshalling.

Type hints will not become mandatory in Python even by convention.

Type declarations are a good place to use the 80/20 rule: Do the easy 20% of the work to find 80% of
the type bugs. Complicated type expressions get challenging to read and write. Mainly, put
type hints on functions and occasionally on variables when the type checker can't figure it out.
E.g. mypy can't infer the element type when you create an empty collection. In a difficult case,
you can punt by using type `Any` which can assign to or from anything.

([These Slides](https://docs.google.com/presentation/d/1xwVHjpQsRTGLdaIpPDNjKPQFtp7mBNfn7oTWK4Kct30/edit?usp=drivesdk)
review our past plan for adding type hints to aid the migration to Python 3.)

You can type a list of strings as `list[str]` or `List[str]`, and similar for other
standard collection types. The capitalized names were created to retrofit type hints into Python 2.
They're unnecessary in Python 3.9+ (wcEcoli is on Python 3.11) but there's no need to change
existing code.


```python
def emphasize(message: str, suffix: str = ''):
  """Construct an emphatic message."""
  # Note that `suffix` is an optional argument to the caller but it's always a
  # `str` to the type checker, so declare it as `str`, not `Optional[str]` (which
  # means `str | None`).
  return message + '!' + suffix
```


## The Typing module

To retrofit static typing analysis to existing releases of Python 2, mypy designers added:
* the type hint comment syntax,
* [the "typing" module](https://docs.python.org/3/library/typing.html) on [PyPI](https://pypi.org)
  containing types like `List` with uppercase names so they didn't have to patch
  existing classes like `list`,
* `.pyi` "stub files" to add type definitions onto existing Python libraries and native libraries,
* the [Typeshed](https://github.com/python/typeshed) repository for type "stub" files.

* `None`: the value None (easier than writing `NoneType`)
* `Any`: any type; assigns to or from anything
* `Dict[str, float]`: a dict of str to float, that is floats keyed with strings
* `Mapping[str, float]`: an abstract mapping; accepts a `dict`, `OrderedDict`, or other
* `List[int]`: a list of int values
* `Tuple[int, int, int]`: a tuple containing 3 integers
* `Optional[str]`: a str or None
* `Callable[[str, List[int]], int]`: a function or anything else you can call given a str and a
list of ints, returning an int
* `str | int`: accepts either a str or an int, aka `Union[str, int]` in code that predates Python 3.10
* `Sequence[float]`: any sequence of floats that supports `.__len__()` and `.__getitem__()`
* `Iterable[float]`: any iterable collection of floats; accepts a `list`, `tuple`, 

* It should work to write numpy type hints like `np.ndarray` and `np.ndarray[int]`.

  ```python
  import numpy as np
  
  def f(a):
      # type: (np.ndarray[float]) -> np.ndarray[int]
      return np.asarray(a, dtype=int)
  ```


### Tools

[PyCharm](https://www.jetbrains.com/help/pycharm/type-hinting-in-product.html) checks types interactively
while you edit. You can also run a particular kind of inspection (such as "Type checker")
or all inspection types on a single file or any source directory. See
[Python Type Checking (Guide)](https://realpython.com/python-type-checking/).

[mypy](https://github.com/python/mypy/) is a batch program to check types. We run it in
Continuous Integration. Code changes must run cleanly to pass the `ecoli-pull-request`
and the `ecoli-small` CI tests. To run it:

    mypy

[PyAnnotate](https://github.com/dropbox/pyannotate) (Py2.7) and
[MonkeyType](https://github.com/Instagram/MonkeyType) (Py3) will observe types at runtime and
write them out as stub files or proposed type hints in the source code.


### Tips

* When a method overrides a superclass method, it inherits the type hints. You needn't
repeat them but if you do, they must match.
* Call `reveal_type(x)` to print type inferences for `x` when mypy runs.
* Escape hatches: `Any`, `cast()`, `# type: ignore`, and `.pyi` stub files (esp. for C extensions).
* Gradually add types to existing code, one file at a time, starting with the most heavily used code.
* Run the type checker in development and in Continuous Integration. Fix its warnings to defend progress.
* We could tell mypy to disallow untyped functions in particular modules.
* The config setting `check_untyped_defs = True` will also check the contents of functions that
  don't have types hints. It might find actual bugs.
* Lambdas don't support type hints, so define a function when you want typing.
* To punt on types for a difficult functions, just leave out the type annotations or add a
  `@no_type_check` decorator. That treats it as having the most general type possible, i.e. each arg
  type is `Any` (except a method's `self` or `cls` arg) and the result is `Any`. `@no_type_check`
  might also disable type inferences and checking inside the function.


### Terminology

* A _type_ is for type checking, variable annotations, and function annotations.
  A _class_ is a runtime thing. Every _class_ acts like a _type_, and there are additional types like
  `str | int` which you can't instantiate like a class.
* Type `A` is a _subtype_ of `B` if it has a subset of the values (classification) **and** a superset of the
methods. A value of a subtype can act like a value from its supertypes.
* The type `From` _is consistent with_ (assignable to a variable of) type `To` if
  * `From` is a subtype of `To`, _or_
  * `From` or `To` is `Any`.
* _Type hints_ and _annotations_ just enable tools: docs, type checkers, IDE code completion, etc.
* Python's type checker is a tool that warns about inconsistent types. It has no impact at runtime.
* _Gradual typing_ means adding type hints to existing code, mixing code
  with and without type hints.
* The type checker can _infer the types_ of local variables, `@property` methods, and more.
* _Nominal typing_ is based on class **names**. Python types are mostly "nominal".
* _Structural typing_ is based on structure such as the duck-type protocol `Typing.Sized`
  which means "any object that has a `.__len__()` method".
* Type information on values is _erased_ at runtime.
* _Covariant types_, _contravariant types_, and _invariant types_: _Variance_ is
  how subtyping between values of a generic type (like `tuple[T, U]`) relates to
  subtyping between its parameter type(s) (`T` and `U`). _Covariant types_ vary
  in the same subtype/supertype direction as their parameters, while
  _contravariant types_ vary in the opposite direction, and _invariant types_
  don't let their parameter types vary at all.


### Covariant and Contravariant types

When variance comes into play, it's puzzling unless you know about it.

See the _covariant/contravariant/invariant type_ definitions, above.
The point is that the rules for generic types need to ensure that assigning a
value to a variable won't violate the variable's declared type.

Let's explain it via examples:

   * ```python
     class Shape: pass
     class Circle(Shape): pass
     s: Shape = Circle()
     c: Circle = Shape()  # mypy type error
     ```

     `Circle` is a subclass of `Shape` and thus a subtype of `Shape`, meaning a
     `Circle` can be substituted for a `Shape`. So you can assign a `Circle` value
     to a `Shape` variable but not vice versa, as far as static type analysis is
     concerned. Similarly, `bool` is a subclass of `int`.

   * ```python
     ts: tuple[Shape] = (c,)
     tc: tuple[Circle] = (s,)  # mypy type error
     ```

     You can assign a `tuple[Circle]` value to an `tuple[Shape]` variable without
     a static type error, but not vice versa.

   * ```python
     s1: Shape = ts[0]
     c1: Circle = ts[0]  # mypy type error
     ```

     Getting an element from a `tuple[Shape]` gets a value that's statically typed
     `Shape`, so it's a static type error to assign this value to a `Circle`
     variable even if it's actually a `Circle` instance.

   * `tuple[T]` is _covariant_ with its element type `T`, i.e. its subtyping
     varies in the same direction as its element type. This works for immutable
     generic types like `tuple` and `frozenset` since their elements
     can only be _sources_; type-compatibile assignment _to_ them doesn't come up.

   * ```python
     def fs(x: Shape) -> str:
        return str(x)
     def fc(x: Circle) -> str:
        return str(x)

     def ffc(f: Callable[[Circle], str]) -> str:
        return f(Circle())

     ffc(fs)  # OK, can pass Circles to fs()

     def ffs(f: Callable[[Shape], str]) -> str:
        return f(Circle())

     ffs(fc)  # mypy type error, static types allow it to pass any Shape to fc
     ```

     `Callable[[Shape], str]` is a subtype of `Callable[[Circle], str]`, meaning
     a function(Shape) is substitutable for a function(Circle), i.e. it's
     OK to assign a function(Shape) value to a function(Circle) variable. Static
     types allow passing Circles to that function -- it's fine receiving
     `Circle` arguments.

   * `Callable[[P1, P2], R]` is _contravariant_ with its parameter types `P1`
     and `P2`, i.e. the subtyping goes the other direction since its parameters
     are _sinks_ of the argument values passed in.

   * ```python
     lc: list[Circle] = [Circle()]
     ls: list[Shape] = lc  # mypy type error
     ls[0] = Shape()
     c2: Circle = lc[0]
     ```

     You might expect a `list[T]` to be _covariant_ with `T`, but there's a problem.
     Such a static type would let us assign a Circle list `lc` to a Shape
     list variable `ls`, then store (or append) any `Shape` into the Shape list
     `ls`. That breaks `lc` -- now the `Circle` list contains a `Shape`, and
     assiging one of those Shapes to a `Circle` variable violates the variable's
     type. To avoid this, it's a type error to assign a `Circle` list to `ls`.

   * Assigning to a `list[T]` variable is _invariant_ with its element type `T`;
     the assigned value must be another list of `T`, no subtypes or supertypes
     of `T`. This rule fits since lists are mutable: they allow getting and
     setting elements.

   * See [Invariance vs
     covariance](https://mypy.readthedocs.io/en/stable/common_issues.html#variance)
     for more info on dealing with type inference (where the type checker infers
     types that weren't declared) and mutable generic collections.

   * `T = TypeVar('T', covariant=True)` declares covariance for a generic type
     w.r.t. its type variable `T`.

   * `T = TypeVar('T', contravariant=True)` declares contravariance.

Also see the short tutorial at [The Ultimate Guide to Python Type
Checking](https://realpython.com/python-type-checking/#covariant-contravariant-and-invariant).

### References

[PEP 484: Type Hints](https://www.python.org/dev/peps/pep-0484/). Worth skimming.

[The "typing" module](https://docs.python.org/3/library/typing.html) defines types like `Any`, `Optional`, `Union`, `Iterable`, `Sized`, and `IO`; and functions like `cast()` and `assert_type()`.

[Python type checking tools](https://python-type-checking.readthedocs.io/en/latest/tools.html)
   * mypy [overview](https://mypy-lang.org/), [docs](https://mypy.readthedocs.org),
   [blog](https://mypy-lang.blogspot.com), [github](https://github.com/python/mypy).
   * PyCharm
     * It can help [add type hints for you](https://www.jetbrains.com/help/pycharm/type-hinting-in-product.html):
     Press Opt-Return (or click the "intention" lightbulb icon) then "Add type hint for ..."
     * [Type Hinting in PyCharm](https://www.jetbrains.com/help/pycharm/type-hinting-in-product.html)
     shows how to use PyCharm features to quickly add type hints, validate them, and eventually convert
     comment-based type hints to type annotations.
   * [pytype](https://github.com/google/pytype) Google's static analyzer checks types, attribute names, and more.
   * [MonkeyType](https://github.com/Instagram/MonkeyType) runtime type info collector


# Performance Tips

[Summarized from sources like [PythonSpeed](https://wiki.python.org/moin/PythonSpeed).]

* Testing membership in a set or a dict is very fast, `O(1)`, unlike a list, tuple, or array.
* Sorting a list using a _sort key_ is faster than using a _comparison function._
* Mapping a function over a list, or using a list comprehension or generator comprehension, should be faster than a `for` loop since it pushes the loop work into compiled C code.
* Local variables are faster to access than global variables, builtins, and attribute lookups.
* Iterators are more memory-friendly and scalable than list operations, so a method like `a_dict.items()` returns an iterator. (Call `list(a_dict.items())` if you need a list.)
* Core building blocks are coded in optimized C, including the builtin datatypes (lists, tuples, sets, and dictionaries) and extension modules like `array`, `itertools`, and `collections.deque`.
* Builtin functions run faster than hand-built equivalents, e.g. `map(operator.add, v1, v2)` is faster than `map(lambda x, y: x+y, v1, v2)`.
* For queue applications using `pop(0)` or `insert(0,v)`, `collections.deque()` offers superior `O(1)` performance over a list because it avoids the `O(n)` step of rebuilding a list for each insertion or deletion.
* Chained comparisons like `x < y < z` are faster and hopefully more readable than `x < y and y < z`.
* Threading can improve the response time in applications that would otherwise waste time waiting for I/O.
