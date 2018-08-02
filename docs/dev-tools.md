Development tools
===================================================

Most of the team uses Sublime Text 3 and PyCharm Pro.

Also install git and maybe the GitHub Desktop app.

On macOS, you might need to install Xcode, then run xcode-select --install (install its command line tools) in order to run the C compiler.

On macOS, iTerm2 is much nicer than the built-in Terminal app.

PyCharm setup
------------------

After building the pyenv and cloning the repo to a local directory, you can create a project in PyCharm.

Select the project's Python interpreter: PyCharm > Preferences > Project: wcEcoli > Project Interpreter > Add > Existing environment > [navigate to something like /usr/local/var/pyenv/versions/wcEcoli2/python].
Set Keyboard Shortcuts: Duplicate one of the available sets (e.g. "Mac OS X 10.5+"), then make changes to suit. A great change is to set Cmd-D (or Ctrl-D) for "Edit > Find > Add Selection for Next Occurrence". [This is like find_under_expand in Sublime Text. Cmd-G (or Ctrl-G) works like Sublime's find_under_expand_skip.]
[TODO] Notes on setting up code style, inspections, ...

Great PyCharm features to know
Cmd-click (or Ctrl-click) a name to jump to its definition.
Navigate > Class, Navigate > File, Navigate > Symbol -- jump to any class, file, or symbol defined in the project or its libraries. This supports "fuzzy matching" where you can type characters in the name to narrow down the list.
Edit > Find > Find Usages -- find all references to the selected symbol.
Double-press Shift -- search everywhere in the project.
Refactor -- commands to do small code refactorings like rename a function or change its calling signature.
Hit F1 to get documentation on an identifier's definition.