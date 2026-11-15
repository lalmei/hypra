"""Generate a categorized GitHub release draft from merged pull requests."""

import argparse
import base64
import json
import re
import subprocess
import tomllib
from urllib.parse import quote, urlencode

START = "<!-- release-drafter:start -->"
END = "<!-- release-drafter:end -->"
VERSION = re.compile(r"\d+\.\d+\.\d+(?:(a|b|rc)\d+)?", re.ASCII)
TITLE_TYPE = re.compile(r"^(?:[^\w\s]\S*\s+)?([a-z]+)(?:\([^)]+\))?!?:\s", re.IGNORECASE)
CATEGORIES = {
    "Features": {"feat", "feature", "enhancement"},
    "Bug fixes": {"fix", "bug", "bugfix", "regression"},
    "Documentation": {"docs", "documentation"},
    "Testing and CI": {"test", "tests", "ci"},
    "Maintenance": {"build", "chore", "dependencies", "maintenance", "perf", "refactor", "style"},
}


def github(endpoint, *, paginate=False, method="GET", payload=None):
    command = ["gh", "api", "--method", method, endpoint]
    if paginate:
        command.extend(["--paginate", "--slurp"])
    if payload is not None:
        command.extend(["--input", "-"])
    result = subprocess.run(
        command,
        input=json.dumps(payload) if payload is not None else None,
        stdout=subprocess.PIPE,
        check=True,
        text=True,
    )
    return json.loads(result.stdout)


def category_of(pr):
    labels = {label["name"].casefold() for label in pr["labels"]}
    for category, names in CATEGORIES.items():
        if labels & names:
            return category

    match = TITLE_TYPE.match(pr["title"])
    if match:
        for category, names in CATEGORIES.items():
            if match[1].casefold() in names:
                return category
    return "Maintenance"


def render_notes(pulls, repo_url, previous_tag, target):
    sections = {category: [] for category in CATEGORIES}
    for pr in sorted(pulls, key=lambda pr: pr["number"]):
        number = pr["number"]
        title = re.sub(r"([\\`*_{}\[\]<>#@&])", r"\\\1", " ".join(pr["title"].splitlines()))
        login = (pr["user"] or {"login": "ghost"})["login"]
        if not re.fullmatch(r"[\w-]+(?:\[bot\])?", login, re.ASCII):
            raise ValueError(f"Invalid author handle for #{number}")
        sections[category_of(pr)].append(
            f"- {title} [#{number}]({repo_url}/pull/{number}) `@{login}`"
        )

    blocks = [
        f"## {category}\n\n" + "\n".join(lines) for category, lines in sections.items() if lines
    ]
    if not blocks:
        blocks.append("No merged pull requests in this range.")
    blocks.append(f"[Compare with {previous_tag}]({repo_url}/compare/{previous_tag}...{target})")
    return "\n\n".join(blocks)


def draft_body(current, notes, version):
    if current is None:
        return f"Hypra {version}\n\n<!-- Add a short feature summary above the PR list. -->\n\n{START}\n{notes}\n{END}\n"
    body = current["body"] or ""
    if body.count(START) != 1 or body.count(END) != 1 or body.index(START) > body.index(END):
        raise ValueError("The existing draft needs one ordered release-drafter marker pair")
    prefix, _ = body.split(START, 1)
    _, suffix = body.split(END, 1)
    return f"{prefix}{START}\n{notes}\n{END}{suffix}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, help="GitHub owner/repository")
    parser.add_argument("--previous-tag", required=True)
    parser.add_argument("--target", required=True, help="Full release commit SHA")
    parser.add_argument(
        "--apply", action="store_true", help="Save the draft instead of printing it"
    )
    args = parser.parse_args()
    if not re.fullmatch(r"[\w.-]+/[\w.-]+", args.repo, re.ASCII):
        parser.error("--repo must be owner/repository")
    if not args.previous_tag.startswith("v") or not VERSION.fullmatch(args.previous_tag[1:]):
        parser.error("--previous-tag must be a Hypra release tag")
    if not re.fullmatch(r"[0-9a-f]{40}", args.target):
        parser.error("--target must be a full lowercase commit SHA")

    repo = f"repos/{args.repo}"
    repository = github(repo)
    default_branch = repository["default_branch"]
    main_comparison = github(
        f"{repo}/compare/{args.target}...{quote(default_branch, safe='')}?per_page=1"
    )
    if main_comparison["status"] not in {"ahead", "identical"}:
        raise ValueError("The target commit must be on the default branch")

    manifest = github(f"{repo}/contents/pyproject.toml?ref={args.target}")
    version = tomllib.loads(base64.b64decode(manifest["content"]).decode())["project"]["version"]
    match = VERSION.fullmatch(version)
    if not match:
        raise ValueError(f"Unsupported Hypra release version: {version}")
    tag = f"v{version}"
    if tag == args.previous_tag:
        raise ValueError("The previous tag must differ from the release version")
    prerelease = match[1] is not None

    matches = [
        release
        for page in github(f"{repo}/releases?per_page=100", paginate=True)
        for release in page
        if release["tag_name"] == tag
    ]
    if len(matches) > 1:
        raise ValueError(f"Multiple releases match {tag}")
    current = matches[0] if matches else None
    if current and (
        not current["draft"]
        or current["target_commitish"] != args.target
        or current["prerelease"] != prerelease
    ):
        raise ValueError(f"{tag} already exists with different release settings or is published")

    previous_sha = github(f"{repo}/commits/{quote('tags/' + args.previous_tag, safe='')}")["sha"]
    pages = github(f"{repo}/compare/{previous_sha}...{args.target}?per_page=100", paginate=True)
    commits = {commit["sha"] for page in pages for commit in page["commits"]}
    if pages[0]["status"] != "ahead" or len(commits) != pages[0]["total_commits"]:
        raise ValueError("Could not verify the complete commit range from the previous tag")

    query = urlencode({"state": "closed", "base": default_branch, "per_page": 100})
    pulls = {
        pr["number"]: pr
        for page in github(f"{repo}/pulls?{query}", paginate=True)
        for pr in page
        if pr["merged_at"] and pr["merge_commit_sha"] in commits
    }
    notes = render_notes(pulls.values(), repository["html_url"], args.previous_tag, args.target)
    body = draft_body(current, notes, version)
    if not args.apply:
        print(body, end="")
        return

    tag_ref = f"refs/tags/{tag}"
    tag_refs = github(f"{repo}/git/matching-refs/tags/{tag}")
    if any(ref["ref"] == tag_ref for ref in tag_refs):
        tagged = github(f"{repo}/commits/{quote('tags/' + tag, safe='')}")
        if tagged["sha"] != args.target:
            raise ValueError(f"Existing tag {tag} points to a different commit")

    payload = {
        "body": body,
        "tag_name": tag,
        "target_commitish": args.target,
        "prerelease": prerelease,
    }
    if current is None:
        saved = github(
            f"{repo}/releases",
            method="POST",
            payload={**payload, "draft": True, "name": f"Hypra {tag}"},
        )
    else:
        endpoint = f"{repo}/releases/{current['id']}"
        latest = github(endpoint)
        fields = ("body", "draft", "tag_name", "target_commitish", "prerelease", "updated_at")
        if any(latest[key] != current[key] for key in fields):
            raise ValueError("The draft changed while generating notes. Rerun to refresh it.")
        saved = (
            latest if latest["body"] == body else github(endpoint, method="PATCH", payload=payload)
        )

    if not saved["draft"] or any(saved[key] != value for key, value in payload.items()):
        raise ValueError("GitHub did not save the expected release draft")
    print(saved["html_url"])


if __name__ == "__main__":
    main()
