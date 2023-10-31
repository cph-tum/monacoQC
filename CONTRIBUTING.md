# Contributing

Feel free to contribute to this project by reporting errors, requesting
features, and of course by fixing and extending the existing source code.

## Report errors or request features

Both can be accomplished by creating an issue on the project's GitHub page.
Please make sure to give a meaningful description of the issue and (in case
of an error) a minimal example that triggers the error. Terminate the issue
title with a dot `.` and describe the necessary steps in order to achieve the
issue's goals.

## Fix or extend existing code

We follow a simplified version of the
[GitLab Flow](https://about.gitlab.com/blog/2023/07/27/gitlab-flow-duo/)
to develop new features and maintain the existing code base. This flow can be
summarized in the following steps:

 - Every non-trivial change starts with an issue. Create one or pick an
   existing one.

 - Create a pull request or merge request with a corresponding feature
   branch. Base the branch on `main` and (if required) adjust its name.

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
      files.", "solver: Fix bug in zero-finding #325." The one-line commit
      message should be limited to 76 characters. Commit messages should
      mention the topic of the commit followed by a colon and an overview
      description with 76 characters. If a detailed description is required,
      add a multi-line commit description. Mention the issue the commit is
      related to with the `#` notation. Every one-line commit message must
      be terminated with a dot `.`.

    - Author name and e-mail are set correctly.

    - Each commit is dedicated to one particular chunk of work.

    - The commits base on the latest main branch so that a fast-forward merge
      is possible.

 - Once you feel your code is ready, assign the merge request to a maintainer
   and remove the `Draft` status of the request.

 - The maintainer may decide to restart the discussion or otherwise performs
   the merge. The feature branch is removed during the merge.

## Formatting

The following rules ensure that redundant white space changes are avoided and
the source code is in a easily readable form:

 - Wrap the code lines in order to improve readability when working in the
   terminal (in particular with git). The general limit is 78 characters.
   However, different languages may use or even require a different limit.
 - Avoid git commits with whitespace changes. The general rule is to have one
   and only one newline at the end of the file and no trailing whitespaces at
   the end of line.
 - Use spaces instead of tabs. Generally, use 4 spaces per indentation.
 - Document your code as you write it. You will not do it later. Do it now.
   Also, do not forget to include pointers to literature (first author, year,
   DOI link) to indicate from which source you took an equation or material
   parameter.
 - Make sure that the files you check in are encoded in UTF-8 and have UNIX
   style line endings (git can take care of the latter, but you may want to
   check if this setting is active).

### Matlab source files

Generally, use small_caps_in_names_with_underscores for names. Exceptions
only where it makes sense, e.g., for the class AlGaAs.

## Platform independence

In order to provide platform independent code, please follow the guidelines
below.

### File name convention

All generated files have lower case file names (including extension). This
facilitates the operation on both case-insensitive and case-sensitive file
systems.

Helper functions such as fullfile(path, file) in MATLAB are recommended as
they automatically provide the correct file path for each operating system.
