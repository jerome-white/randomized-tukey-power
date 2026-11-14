# Contributing

This file records the conventions that contributors, human or
otherwise, are expected to follow. It is short on purpose: where a
convention is the standard one in the Julia ecosystem it is named, not
re-explained. It is adapted from the conventions of
[RandomizedTukeyHSD.jl](https://github.com/jerome-white/RandomizedTukeyHSD.jl),
whose power this repository studies.

## Making a change

These rules apply to every contributor, human or agent:

- **Every change has an issue**; however small. Describe the problem
  and the intended change fully enough that someone without your
  context can act on it. Length should be commensurate with achieving
  that level of clarity.
- **One issue, one branch, one pull request.** Branch from `master`
  and name the branch `<type>/<issue>/<short-slug>`. This extends
  [Conventional Branches](https://conventionalbranch.org/): `<type>`
  is one of `feature`, `bugfix`, `hotfix`, `release` or `chore`, and
  `<short-slug>` is a few lowercase words joined by hyphens. We add
  `<issue>`, the number of the issue, as its own segment, which the
  convention does not define. Examples: `feature/12/step-down-method`
  or `bugfix/30/seed-reuse`. Use `chore` for changes to
  documentation, tooling and other upkeep. The pull request description
  starts with `Closes #N`, explains what changed, why it changed, and the
  problem it solves. It should also say what was verified and what was
  not.
- **Commit at each logical step**; not once at the end. A commit
  message has a summary line, a body that says what changed and why,
  and ends with `Refs #N` (where N is the issue being
  addressed). Commit history is viewed as part of the project
  documentation: `git blame` should answer "why is it like this".
- **The maintainer merges.** Do not push to `master`, and do not merge your
  own pull request.
- **Issues are grouped into families.** A larger thread of work is a parent
  issue, and the issues that carry it out are its sub-issues. Open a
  sub-issue under the parent it belongs to, and say in its text how it
  contributes to the parent.

## Style

No formatter is enforced, so please stick to the following rules manually:

- Lines are at most 80 columns, and code is ASCII.
- Check an argument with an `if` block that throws, not with
  `cond || throw(...)`.
- No module-level mutable globals.
- A helper used once is written inline, not as a function. The one
  exception is a script's `main`, which keeps its loops out of global scope.

In general, follow the style of the file you are editing. The code is
Julia.

## Research conventions

This repository is a research project, so some conventions concern
evidence rather than code.

- **Every claim has a status.** *Verified* means checked against a
  source that was read, or against a run whose command and output are
  recorded, as described under "Every number that is cited can be
  re-derived".
  *Reasoned* means it follows from an argument but has not been checked.
  *Unverified* means it comes from memory, conversation or hearsay. A claim
  moves up a level only with evidence, and the evidence is cited where the
  claim is made. When in doubt, use the lower level.
- **Experiments are reproducible.** A run is determined by a seed, its
  parameters and the commit of the code. A result states all three, so
  that anyone can regenerate it. The overarching principle is that the
  process of obtaining a result is as important as the result itself.
- **Every number that is cited can be re-derived.** The test is whether,
  if this session disappeared, someone could get the number again from
  what is in the repository. Exploration is encouraged: an agent may try
  things as the maintainer would. What is quoted in an issue, pull request
  or document, or used to choose a direction, is recorded where it is
  quoted: the command, the seed, the commit and the output. A script is
  committed when it will be rerun or belongs to the experiment suite; a
  one-off check leaves its command in the comment. A value that is only
  estimated, such as one read off a plot, is labelled as an estimate.
- **Fix the plan before looking at results.** The question, the metrics and
  the rule for deciding what counts as a win are written in an issue before
  the experiment that tests them runs. If the plan changes after results
  are seen, the change and its reason are recorded, not quietly applied.
- **Long runs belong to the maintainer.** A job that takes more than a
  minute or two is run by the maintainer. An agent's commands have time
  limits, and a run the maintainer performs leaves an explicit, repeatable
  instance of how the result was made. The agent supplies the command, the
  reason for it, the expected duration and the output path. An agent starts
  no background jobs.
- **Discuss the science.** Issues, pull requests, commits and documents are
  about the questions, methods and evidence. Where the work will be
  published is not discussed here.

## Use of AI

This repository is developed with an AI coding agent, and contributors may
use one too. This section says how, so that readers, reviewers and journals
can judge the work, and what is expected of anyone who contributes with AI.

### How this repository is developed

- Work began on 2026-11-13 with Claude Code, using Claude Sonnet 5.5
  (`claude-sonnet-5-5`). The maintainer sets the direction and the
  research questions. The agent drafts the issues, documents, code, commit
  messages and pull requests. The maintainer reads them on GitHub, asks for
  changes, and merges every pull request.
- Every commit an agent wrote has a `Co-Authored-By:` trailer naming the
  model, so the model behind any commit is in the history:

  ```bash
  git log --format='%h %(trailers:key=Co-Authored-By,valueonly)'
  ```

  Commits that do not contain that trailer came from a human.
- Review is not the only evidence of correctness. A result is trusted to the
  degree that it is reproducible from its seed and recorded parameters, and
  that each claim carries the status described above.

### What is expected of contributors

These are requirements. They apply to both people and agents.

- **Say that AI helped**, in the pull request description. Keep the
  `Co-Authored-By:` trailer on commits an agent wrote.
- **You are accountable for what you post** (issues, comments, commits,
  pull requests, code), whoever or whatever wrote it. Maintainers may close
  contributions with unverified claims presented as verified, text that
  does not match the code or the data, or changes nobody asked for.
- **Do not run an agent unattended.** An agent works in a session that a
  person is running and answers for.
- **Hold AI-written code and documentation to the same standard.** A test
  must test something, and a description must be true of what it describes.
- **Say what you verified** in the pull request, and what you did not.
- The rules in "Making a change" apply, including an issue for every change.
