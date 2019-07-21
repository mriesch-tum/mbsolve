# Contributing

Feel free to contribute to this project by reporting errors, requesting
features, and of course by fixing and extending the existing source code.

## Report errors or request features

Both can be accomplished by creating an issue on the project's GitHub or
GitLab page, respectively. Please make sure to give a meaningful description
of the issue and (in case of an error) a minimal example that triggers the
error.

## Fix or extend existing code

We follow the
[gitlab flow](https://docs.gitlab.com/ee/workflow/gitlab_flow.html)
to develop new features and maintain the existing code base. This flow can be
summarized in the following steps:

 - Every non-trivial change starts with an issue. Create one or pick an
   existing one.

 - Create a pull request or merge request with a corresponding feature branch.
   Base the branch on "master" and (if required) adjust its name.

 - Work on the solution to the issue and push to the branch to share and
   discuss your code changes with the other developers.

 - Feel free to rearrange the commits on the feature branch (e.g., rebase
   with squash). In fact, it is your duty to provide clean and non-redundant
   commits before your code can be merged.

   Please note that it is not allowed to alter commits outside of the feature
   branches. Also, be careful when (or avoid altogether) basing your work
   on a feature branch.

   Clean commits:

    - Meaningful commit message, e.g., ".gitignore: Added generated output
      files.", "monte-carlo: Fixed indentation." The one-line commit message
      should be limited to 76 characters. If a detailed description is
      required, add a multi-line commit description.

    - Author name and e-mail are set correctly.

    - Each commit is dedicated to one particular chunk of work.

    - After each commit all continuous integration (CI) jobs pass, in
      particular the "format" job (see "Formatting" below).

    - The commits base on the latest "master" branch so that a fast-forward
      merge is possible.

 - Once you feel your code is ready, assign the request to a maintainer and
   remove the WIP status of the request.

 - The maintainer may decide to restart the discussion or otherwise performs
   the merge. The feature branch is removed during the merge.

## Formatting

Please follow the rules below when contributing to the mbsolve project.

### General rules for all source code and text files

 - Avoid lines with more than 78 characters.

 - Remove whitespaces at the end of the line. However, insert a newline at
   the end of the file.

 - Use UNIX-style newlines (In order to use e.g. Windows-style newlines
   locally, you can adjust git accordingly).

 - Use snake case for names (lower case and separate words with underscores,
   e.g. my_class).

 - Keep names as short as possible, yet provide a meaningful name (do not
   name your variables "counter_of_things_i_count", neither name it "ctrtic"
   -- consider refactoring your code so that you can use "counter" without
   any collision).

### Special rules for CMake source code files

 - Use two spaces for indentation. No tabs.

 - Follow the conventions used in the existing code.

### Special rules for C++ source code files

 - The code formatting is based on the [Mozilla Style Guide](https://developer.mozilla.org/en-US/docs/Mozilla/Developer_guide/Coding_Style)
   The points where this style is overruled are defined in this file.

 - Use four spaces for indentation. No tabs.

 - Use C-style comments even for single lines (/* like this */)

 - Group the include statements:
    1. Includes of standard library components
    2. Includes of system libraries, compilers, external frameworks
    3. Includes of the Eigen library
    4. mbsolve headers

   and sort each group. Use the version with angular brackets.

       #include <solver.hpp> /* Good */
       #include "solver.hpp" /* Bad */

 - The constructor initializers should look like

       my_class(many, arg, here)
         : m_many(many), m_arg(arg),
           m_here(here)

 - Do not indent case labels.

 - After C-style casts, there is a space.

       (int) i /* Good */
       (int)i  /* Bad */

 - In case the arguments do not fit into a single line, arrange them such
   that every argument gets a single line:

       function_call(
           arg1,
           arg2,
           another_function_call(arg3),
           arg4);

 - Do not use nested namespaces, use only the namespace "mbsolve".

 - Run clang-format to format the code automatically. It makes sense to
   integrate this tool into your editor. At least before issuing a
   pull request it is mandatory that the code is in agreement with the
   rules above.

### Special rules for Python source code files

 - Use four spaces for indentation. No tabs.

 - Follow the conventions used in the existing code.

### Tools

#### editorconfig

The file .editorconfig is provided for editors that use the editorconfig
project.

#### EMACS

Add the following to your .emacs file:

    ;; indent with spaces
    (setq-default indent-tabs-mode nil)

    ;; delete trailing whitespace when saving
    (add-hook 'before-save-hook 'delete-trailing-whitespace)

    ;; maintain newline at the end of the file
    (setq require-final-newline t)

    ;; treat CUDA code as C++ source files
    (add-to-list 'auto-mode-alist '("\\.cu$" . c++-mode))

    ;; treat template files as C++ source files
    (add-to-list 'auto-mode-alist '("\\.tpp$" . c++-mode))

    ;; indentation size 4 spaces
    (setq-default c-basic-offset 4)

    ;; indent after open brace for parameters
    (c-set-offset 'arglist-intro '+)

    ;; do not indent within namespace region
    (c-set-offset 'innamespace 0)

    ;; automatically reload buffers
    (global-auto-revert-mode t)

#### clang-format

From the command line, run

    $ clang-format -i file.cpp

to format file.cpp in place according to the rules given in the .clang-format
file.

clang-format can be integrated easily into several IDEs and editors. For
example, in EMACS it can be executed at save

    ;; run clang-format when saving in C++ mode
    (load "[path-to-clang-format.el]/clang-format.el")
    (add-hook 'c++-mode-hook
              (lambda ()
                (global-set-key [f8] 'clang-format-buffer)
                (add-hook 'before-save-hook 'clang-format-buffer nil 'local)))

## Documentation

 - Document your code using Doxygen annotation.
