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

## 2. Deploy it (30 minutes, optional but worth a prize)

There is a separate prize — "Best Use of Deployxa" — and the project page is
much stronger with a URL a judge can click than with a `git clone` instruction.

The container needs no configuration: no API keys, no database, and the ClinVar
baseline ships inside the image.

```bash
docker build -t dissent . && docker run -p 8000:8000 dissent
```

That is verified to work, including under a read-only filesystem. For Deployxa
or any similar platform, point it at the repository; it will find the
`Dockerfile`. The only setting worth adding is a writable volume for the cache:

```bash
DISSENT_CACHE_PATH=/var/cache/dissent/dissent.sqlite
```

Without it the server still runs — the cache turns itself off and
`/api/health` says so — it is just slower.

After deploying, check these four by hand: `/`, `/panel`, `/scale`, and one
real variant typed into the box.

Then tick "Best Use of Deployxa" on the Devpost form.

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
- 519 tests, passing with the network severed at the socket layer; CI runs them
  on Python 3.11, 3.12 and 3.13
- Verified from a fresh `git clone` into an empty directory: install, tests,
  CLI, panel CLI, and the web server all work with no configuration
- Dockerfile built and run, including under `--read-only`, with graceful
  shutdown confirmed
- Every page meets WCAG AA contrast in both themes, with keyboard focus rings
  and skip links
- All Devpost prose written, and the measured numbers in it checked against the
  code by the test suite
