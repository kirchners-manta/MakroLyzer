# Issue tracker: Local Markdown

Issues and specs live in `.scratch/`, which is temporary, gitignored
working space.

- One feature per directory: `.scratch/<feature-slug>/`.
- Spec: `spec.md` within that directory.
- Tickets: `issues/<NN>-<slug>.md`, numbered from `01`,
  with one file per ticket.
- Record triage state near the top as `Status: <value>`,
  using `triage-labels.md`.
- Append discussion under `## Comments`.

When a skill says to publish to the issue tracker, create the
appropriate local file. When it says to fetch a ticket, read the
referenced file. Resolve ambiguous ticket numbers within their
feature directory.

Create directories and files only when there is work to record.
