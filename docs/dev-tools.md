# Development tools

Most of the team uses **Sublime Text 3** and **PyCharm Pro**.

Also install **git** and maybe the **GitHub Desktop app**.

On macOS, you need to have Xcode's command line tools to run the C compiler. Unless homebrew did it for you, run

    xcode-select --install

**Tip:** On macOS, **iTerm2** is much nicer than the built-in Terminal app.

## Additional steps for macOS Mojave

**Mojave:** **Before upgrading to macOS Mojave, update FUSE.** After installing Mojave, run `xcode-select --install` again.

**Mojave:** For macOS Mojave or higher (10.14+), follow the "you will also need to install the additional SDK headers" instructions on [https://github.com/pyenv/pyenv/wiki/Common-build-problems](https://github.com/pyenv/pyenv/wiki/Common-build-problems). In short:

   ```bash
   sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /
   ```

## Installing the pyenv tools

pyenv lets you easily switch between multiple versions of Python and between groups of third party libraries.

1. Install `pyenv`, `pyenv-virtualenv`, `pyenv-virtualenvwrapper` using your local package manager, e.g. **homebrew** on macOS. (These are already installed on Sherlock.) E.g.
   ```bash
   brew install pyenv pyenv-virtualenv pyenv-virtualenvwrapper
   ```
2. Initialize `pyenv` and optionally `pyenv-virtualenv` in your shell login script (`.bash_profile` on Linux; `.profile` or `.bash_profile` on macOS, etc.).
   - See [setup_using_python.md](https://github.com/CovertLab/ComputationalResources/blob/master/_sherlock/setup_using_python.md) and docs in the same directory about setting up this stuff including Linux modules, but beware that those docs are out of date, e.g. they're for now-obsolete Sherlock 1.0.

   - Example `.profile` or `.bash_profile` lines for macOS:

   ```bash
   export PYENV_ROOT=/usr/local/var/pyenv
   if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi
   if which pyenv-virtualenv-init > /dev/null; then eval "$(pyenv virtualenv-init -)"; fi
   ## -- Do this before sourcing iterm2_shell_integration
   ```

   - Example `.bash_profile` lines for Sherlock:

   ```bash
   module load wcEcoli/sherlock2

   export PYENV_ROOT="${PI_HOME}/pyenv"

   if [ -d "${PYENV_ROOT}" ]; then
       export PATH="${PYENV_ROOT}/bin:${PATH}"
       eval "$(pyenv init -)"
       eval "$(pyenv virtualenv-init -)"
   fi
   ```

3. You'll need to put the project on the `PYTHONPATH` when working on it. Consider adding this to your profile or creating a shell alias to do this when you work on wcEcoli:

   ```bash
   export PYTHONPATH="$HOME/wcEcoli:$PYTHONPATH"
   ```

4. Open a new shell so it runs the updated profile.

5. NOTE: if you have a `~/.local/` directory, paths might not work properly with `pyenv` and you might receive error messages. ([TODO] In that case, delete the directory?)


## PyCharm setup

After building the pyenv and cloning the repo to a local directory, you can create a project in PyCharm.

* Select the project's Python interpreter: PyCharm > Preferences > Project: wcEcoli > Project Interpreter > Add > **Existing environment** > [navigate to something like `/usr/local/var/pyenv/versions/wcEcoli2/python`].
* Set Keyboard Shortcuts: Duplicate one of the available sets (e.g. "Mac OS X 10.5+"), then make changes to suit. A great change is to set Cmd-D (or Ctrl-D) for "Edit > Find > Add Selection for Next Occurrence". [This is like `find_under_expand` in Sublime Text. Cmd-G (or Ctrl-G) works like Sublime's `find_under_expand_skip`.]

[TODO] Notes on setting up code style, inspections, ...

### Great PyCharm features to know

* Cmd-click (or Ctrl-click) a name to jump to its definition.
* Navigate > Class, Navigate > File, Navigate > Symbol -- jump to any class, file, or symbol defined in the project or its libraries. This supports "fuzzy matching" where you can type characters in the name to narrow down the list.
* Edit > Find > Find Usages -- find all references to the selected symbol.
* Double-press Shift -- search everywhere in the project.
* Refactor -- commands to do small code refactorings like rename a function or change its calling signature.
* Hit `F1` to get documentation on an identifier's definition.


## iTerm2 Tips

**iTerm2** is a macOS app that [greatly improves on the stock terminal app](https://www.iterm2.com/features.html).

Tips (also see [the FAQ](https://www.iterm2.com/faq.html)):

* If you configure it to save & load preferences to a Dropbox folder, you don't have to do much when switching to a new Mac.
* tmux Integration lets you make and adjust window split panes much more easily than typing tmux commands.
* [Shell Integration](https://www.iterm2.com/documentation-shell-integration.html) is very handy, but the regular setup didn't work quite right on Sherlock with the pyenv virtualenv shell prompt. So for Sherlock, Jerry just set up the Triggers as documented on that page. The "Prompt Detected" trigger is probably the most useful part since it lets you jump between shell prompts in the terminal output.

Example "Default" profile configuration for "Keys":
* Option up/down arrows: scroll one line up/down
* Shift left/right arrows: send hex code `0x02` or `0x06` (respectively) to move the cursor
* Control left/right arrows: Send `^[b` or `^[f` (respectively; `^[` is ESC) to move the cursor by words (assuming your `.inputrc` file uses EMACS style readline editing, which is the default)
* Option left/right arrows: Same as Control left/right arrows
* Command left/right arrows: send hex code `0x01` or `0x05` (respectively) to move the cursor to the start/end of line
* Left Option (alt) key: Esc+
* Right Option (alt) key: Normal

Example "Default" profile configuration for "Triggers":
* `@login.sherlock.stanford.edu's password:` --> Open Password Manager
* `Enter passphrase for key '/home/users/\w+/.ssh/id_rsa':` --> Open Password Manager
* `^\[(\w+)@(sh-[\w-]+) login! ([^\]]+)]\$ ` --> Prompt Detected
* ^ You can use the same regex for Report User & Host `\1@\2` and for Report Directory `\3`

Example overall configuration for "Keys":
* Control Tab, Control-Shift Tab: Next Tab, Previous Tab
* Shift up/down arrows: Scroll One Line Up/Down
* Command up/down arrow: Scroll to Top/End
* Control up/down arrows, Page Up/Down, Command Page Up/Down, Shift Page Up/Down: Scroll One Page Up/Down
* Command Home/End: Scroll to Top/End
