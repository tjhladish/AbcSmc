# Style

## General Naming

 - Use indicative variable names.
 - Use `pot_hole_case` for function names, arguments, internal variables, etc
 - for variables / functions / etc that are private or not exported, start them
 with an underscore. E.g. `_some_internal_variable` or `_unexported_function_def`
 - There are some shorthand conventions:
  * parameter: `par`
  * metric: `met`
  * index: `idx`

## Global Constants, Macros

These should be UPPER_CASE.

## Structural

 - Prefer use `auto` where practical. While there are advantages to an
 obvious type declaration, modern IDEs can infer type with static compilation
 tools and provide it via mouseover. Using `auto` generally makes refactoring
 smoother.
 - To the extent practical, keep declarations within the `ABC::` namespace.
 - While there are advantages to header-only libraries, AbcSmc practically
 cannot be header-only. Therefore, headers should minimize implementation
 details to ensure minimal implementation leakage, maximum portability, best
 compile times, etc.

## Concept

In general, an `AbcSmc` object is intended to be what end users construct and
interact with. In turn, that object will construct and interact with other
objects that take care of the various tasks in an analysis, according to what
the user supplies in a configuration file.

To the extent possible, separate concerns have their own objects. At one time,
for example, the `AbcSmc` object had all of the details of implementing database
storage and interaction. Now, `AbcSmc` has a storage object that it interacts
with in a generic way. The storage object is responsible for the details of e.g.
interacting with a database _or any other storage mechanism_.