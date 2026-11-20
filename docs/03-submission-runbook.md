# Submission runbook

The parts that need your accounts, in the order to do them. Everything that
could be done without you is done; the text to paste is in
[02-submission.md](02-submission.md).

Official close: **26 December 2026, 00:00 Eastern Time**, which is the *end of
25 December*. Aim for 24 December.

---

## 1. Push the code (15 minutes)

A judge needs a public repository link. The repository is already a clean git
history with no secrets and no large binaries: 2.4 MB of tracked files, of
which 1.3 MB is the measured ClinVar baseline that has to ship.

Create an empty repository on GitHub named `dissent` — **public**, with no
README, no .gitignore and no licence, because the repository already has all
three and GitHub's versions would collide.

Then, from the project folder:

```bash
git remote add origin https://github.com/<your-username>/dissent.git
```

```bash
git push -u origin main
```

If the push is rejected because the branch is called something else, check with
`git branch --show-current` and push that name instead.

**Before pushing, know what becomes public.** Git built the commit author
address from this machine's hostname, so every commit carries it. Check what it
is with `git log -1 --format='%an <%ae>'`. Changing it means `git config
user.email` plus a rewrite of the existing history, so if you want it changed,
do it *before* the first push rather than after. (It is printed here by a
command rather than written out, so that this file does not publish it a
second time.)

Then verify what a judge will actually see:

```bash
gh repo view --web
```

## 2. Deploy it — optional (20 minutes)

**This is not a submission requirement.** The rules list five things to
submit: a project page, the video, a code link, a track, and the team members.
A live URL is not among them, and the only "does it run" requirement is the
video's — *"Show the project working."*

Deploying buys exactly two things:

1. **The "Best Use of Deployxa" prize.** It is specific to that platform —
   *"Opt in when you submit if you deployed your project on Deployxa"* — so
   deploying somewhere else does not qualify.
2. A judge being able to click and try it, which may help Presentation and
   Technical execution. May, not must.

It costs one free account. There is no way around that: every host needs an
account, and the one that is already signed in here, Google Cloud, needs
billing enabled on it before Cloud Run will accept anything — a bigger
commitment than a free signup, not a smaller one.

If you skip this, nothing else in the submission changes. Delete this section
from your plan and go to the video.

Deployxa's free plan fits this project with room to spare: **no credit card,
no time limit, 3 projects, 512 MB RAM** and a free `*.deployxa.online`
address.

### Rehearsed, so you are not debugging on their servers

Everything below except the account has been run end to end: the public
repository was cloned the way the platform clones it, built from the
Dockerfile alone, and run capped at **512 MB and one CPU** — the free plan's
limits.

| Check | Result |
|---|---|
| Builds from the Dockerfile with no configuration | yes |
| All routes serve (`/`, `/panel`, `/scale`, the APIs) | 200 |
| A live variant query end to end | HFE C282Y → Benign, 5 conflicts |
| Peak memory under load (panel + 40 concurrent page loads) | **66 MB of 512** |
| Killed for memory at any point | no |

So the free tier is not a squeeze; it is 13% used.

### If you do deploy: the three steps that need you

**1.** Create the account at <https://deployxa.com> — "Deploy free", no card.

**2.** Dashboard → **New project** → **Git URL**, and paste:

```text
https://github.com/candicewu0515/dissent.git
```

Branch `main`. Use **Git URL** rather than the CLI: the CLI uploads the whole
working folder and has no ignore file, which would send the virtual
environment and the local cache. The repository is public, so no GitHub
connection is needed either.

**3.** Deploy. Deployxa detects the Dockerfile and builds from it — their own
documentation says *"Projects with a Dockerfile build from it and skip those
fields"* — so leave the build command and output directory alone.

**There are no environment variables to set.** The image already points the
cache at a writable volume it owns, and the app reads `PORT` from the
environment. If the platform gives it a read-only filesystem anyway, the cache
turns itself off and the server keeps serving; `/api/health` reports it.

### One platform quirk, found by testing

Deployxa's edge returns **403 to requests whose User-Agent is
`Python-urllib/…`**, the Python standard library's default. Everything else
gets through: curl, `requests`, browsers, a custom agent string, even no
User-Agent header at all. It is a bot filter on their side, not a fault in the
app — the same request succeeds locally.

It does not affect the submission, and no judge will hit it. It is recorded
here only so that a script that suddenly 403s is not mistaken for the site
being down.

### After it is live

Check these four by hand, because a build that succeeds is not the same as an
app that works:

- [ ] `/` — type `HFE C282Y` and get a report
- [ ] `/panel` — load the carrier panel and let it finish
- [ ] `/scale` — the four figures and the gene chart render
- [ ] `/api/health` — says `"status":"ok"`

Then tick **"Best Use of Deployxa"** on the Devpost form, and put the URL in
the project page's "Try it out" links.

## 3. Record the video (half a day, including retakes)

This is 25% of the score on its own, equal to the code. The shot-by-shot script
is at the end of [02-submission.md](02-submission.md).

Before you hit record:

```bash
python analysis/warm_cache.py
```

It fetches everything the demo touches and exits non-zero if any demo query
would fail live, so it can gate the recording. Entries last a week.

Hard requirements from the rules, each of which disqualifies the video on its
own:

- **3 minutes maximum.**
- **YouTube or Vimeo** — not a file upload.
- **Embedding turned on.** This is off by default on some YouTube accounts.
  Set it before you submit and then open the Devpost preview to confirm the
  player actually appears.
- It has to **show the project working**, not describe it.

Captions throughout: several judges will watch muted.

## 4. Fill in the Devpost page (45 minutes)

Paste from [02-submission.md](02-submission.md). The fields map like this:

| Devpost field | Section to paste |
|---|---|
| Project name | "Project name" |
| Elevator pitch | "Elevator pitch" |
| About the project | "What it does", "Why it matters", "How we built it", "Challenges we ran into", "What we learned", "What's next" |
| Built with | the tools list under "Required disclosures" |
| Try it out links | the live URL and the repository |

Two things the rules require that are easy to miss:

- **List the AI tools used.** "You can use open-source libraries, APIs,
  templates and AI tools **if you list them in your project description**."
  The disclosure paragraph is written; paste it as-is.
- **Add every team member on Devpost.** If you are solo, there is nothing to
  do, but confirm the submission shows one participant rather than zero.

Track: **Healthcare & Environmental**. The reasoning is in
[02-submission.md](02-submission.md#track).

## 5. Final check before you submit

Open your own project page as a stranger would, in a private window, and
confirm:

- [ ] The video plays, inline, and is under 3:00
- [ ] The repository link opens without a login
- [ ] The live URL loads and one variant query works
- [ ] The track is set
- [ ] The AI-tools disclosure is visible in the description
- [ ] Submitted, not just saved as a draft

The last one is the single most common way a finished hackathon project scores
zero.

---

## What is already done

Nothing below needs you.

- A clean git history, no secrets, nothing uncommitted
- 537 tests, passing with the network severed at the socket layer; CI runs them
  on Python 3.11, 3.12 and 3.13
- Verified from a fresh `git clone` into an empty directory: install, tests,
  CLI, panel CLI, and the web server all work with no configuration
- Dockerfile built and run, including under `--read-only`, with graceful
  shutdown confirmed
- Every page meets WCAG AA contrast in both themes, with keyboard focus rings
  and skip links
- All Devpost prose written, and the measured numbers in it checked against the
  code by the test suite
