# AI4Transcriptomics 2026

Working repository for the practical sessions of *From Nucleus to Neural
Networks: Foundations of AI in Transcriptomics*, AI4Health@DAISY group, LIG,
Université Grenoble Alpes.

In the practical you redo parts of a published analysis of FSHD muscle biopsies
(Wang et al., Hum Mol Genet 2019, GEO series GSE115650) from the counts the
authors deposited. This page tells you how to set up, how to submit and how you
are graded.

| Where | What you find there |
|-------|---------------------|
| [practical/README.md](practical/README.md) | the four tasks, the dataset, the paper, the scripts and how to run them |
| [practical/worksheets/](practical/worksheets) | one form per step: fill them as you work, your report cites them |
| UNESS | the report template, and where you submit your report |
| [Issues](https://github.com/AI4Health-DAISY-LIG/AI4Transcriptomics_AI4OH_2026/issues) | questions and answers, opened with the question template |

# Practical sessions: instructions and grading

Course repository: https://github.com/AI4Health-DAISY-LIG/AI4Transcriptomics_AI4OH_2026

Final deadline: **Wednesday 16 December 2026, 11:59pm (Paris time zone)**

## 1. How the practicals are graded

- Report and code: 85 % of the practical grade.
- Collaboration: 15 % of the practical grade, detailed in section 5.
- Bonus: up to 2 points for an optional part on Biolink and Translator,
  detailed in section 6.

## 2. Set up (before the first practical, as shown in course 3)

1. Create a GitHub account if you don't have one, and send me your GitHub handle
   through the set-up form on UNESS. I will create a branch named after
   your handle.
2. Fork the course repository to your own account (the "Fork" button, top right
   of the repository page).
3. Clone your fork to your computer and work there (allow for contributions to
   the main repo). Commit and push regularly: your commit history is part of
   your work. Commit often: it is fine if your code is broken or has mistakes in
   it. Each time you end a coding session, commit to your repo with a useful
   message (for example "corrected ...", "to do next time: ...", "having a bug
   ...").

The course repository is public, in line with open science practice: your fork,
code and comments will be visible to anyone. You may use a GitHub account that
does not show your real name. Your handle is only linked to your name on UNESS.

## 3. Submit your work

1. On GitHub, open a pull request from your fork to the course repository.
2. When opening it, check the target: the base repository must be the course
   repository, and the base must be the branch named after your handle, not
   `main`.
3. Open the pull request early. Everything you push to your fork afterwards
   appears in it automatically. The state of your pull request at the deadline
   is what will be graded; commits pushed after the deadline will not be taken
   into account.
4. Your pull request must contain all your code, your results and your filled
   worksheets. It must not contain the data,
   the paper or your report: git ignores the `data/` and `doc/` folders of the
   practical for that reason.
5. Submit your report as a PDF on UNESS before the same deadline: three pages,
   written from the template you find on UNESS. The
   report does not go into the repository. Each factual claim in it must link to
   your code, your results or your worksheets: see
   [practical/README.md](practical/README.md#the-report).

To ask a question, open an Issue on the course repository. The question
template is the only way to open one.

## 4. Collaborate

- Commit to helping 2 or 3 peers, and send me their names or GitHub handles
  through the set-up form. I will make sure every student has at least 2
  committed helpers. Committing means answering when they ask, not doing their
  work.
- Ask and answer questions in the Issues of the course repository, not in
  private messages. One issue per question; link to the relevant file or line of
  code. This is your help trail, and it is what makes your help visible.
- Review at least one peer's pull request in the week before the final deadline:
  comment on their code, suggest improvements, ask questions. Use the feedback
  you receive to finalize your own report.

Allowed and encouraged:

- discussing the paper and the methods,
- explaining concepts,
- debugging together,
- reviewing each other's code,
- pointing to documentation or examples.

Not allowed:

- copying code or text from another student,
- letting someone copy yours (send me a private email to report it).

Substantially identical submissions will be treated as plagiarism. If someone
asks you to do their work, you may refuse and tell me privately. This will not
count against you.

## 5. How collaboration is graded (15 % of the practical grade)

**Milestones (half of the collaboration grade).** Each milestone is checked on
GitHub and is either met or not:

| By | Milestone |
|----|-----------|
| Wednesday 25 November | fork created, first commit pushed, pull request opened to your branch |
| Wednesday 2 December | at least one Issue opened or answered on the course repository |
| Friday 4 December | a first complete draft (code) in your pull request, so that peers can review it |
| Wednesday 9 December | at least one review of a peer's pull request |
| Wednesday 16 December | peer ratings submitted through the peer review form on UNESS |

**Help given and peer feedback (half of the collaboration grade).** With your
final submission, you will privately rate each peer you actually worked with on
two questions, from 1 to 4:

- Did this person help you when you asked? (1 = never, despite being asked;
  4 = always, and the help was useful)
- Did this person contribute their own work? (1 = mostly relied on others;
  4 = clearly did their own work and shared it)

Ratings are averaged: being rated by many peers does not raise your grade, being
genuinely helpful does. I will cross-check them against the issues, pull-request
reviews and commit histories; ratings not supported by this trail will not be
taken into account. Students whose trail shows substantial help to others will
be credited for it even if they received few ratings.

You will not be penalized if a peer you committed to never needed help. Doing
nothing, or relying on others without contributing yourself, will lower this
grade.

## 6. Bonus (up to 2 points)

This part is optional. You take genes from your own results and query the NCATS
Biomedical Data Translator about them, then read one answer as a Biolink
statement.

- Use the production instance only: https://ui.transltr.io/
- Log in with your GitHub account.
- Save each query and keep the link from "Share Result Set".
- Report your queries and their results in the bonus part of your report, on
  UNESS.

The steps are in
[practical/README.md](practical/README.md#bonus-biolink-and-translator).

## License

The code in this repository is under the MIT license (see [LICENSE](https://github.com/AI4Health-DAISY-LIG/AI4Transcriptomics_AI4OH_2026/blob/main/LICENSE)).
# Init TP
