"""Inspect generated docs; browser rendering remains a separate smoke test."""
import re
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / 'site'


def main():
    required = [
        'index.html', 'python-api/index.html', 'cli-reference/index.html',
        'latex-and-docs/index.html', 'architecture/index.html',
        'assets/katex/katex.min.js', 'assets/katex/katex.min.css',
        'assets/mermaid/mermaid.esm.min.mjs', 'mathematics.pdf',
    ]
    missing = [name for name in required if not (SITE/name).is_file()]
    if missing:
        raise SystemExit(f'Missing generated files: {missing}')
    api = (SITE/'python-api/index.html').read_text()
    cli = (SITE/'cli-reference/index.html').read_text()
    if 'hypra.api.expand_epsilon' not in api or 'doc-signature' not in api:
        raise SystemExit('mkdocstrings API reference did not render')
    if 'backend-info' not in cli or '--latex' not in cli or '::: mkdocs-typer2' in cli:
        raise SystemExit('Real Typer command documentation did not render')
    for path in (SITE/'assets/mermaid').glob('*'):
        if path.suffix not in ('.js','.mjs'):
            continue
        for dependency in re.findall(r'(?:from\s*|import\s*(?:\(\s*)?)["\']\./([^"\']+\.(?:m?js))["\']',path.read_text()):
            if not (path.parent/dependency).is_file():
                raise SystemExit(f'Missing offline Mermaid dependency: {dependency}')
    class Assets(HTMLParser):
        def handle_starttag(self, tag, attrs):
            attributes = dict(attrs)
            url = attributes.get('src', '') if tag == 'script' else ''
            if tag == 'link' and attributes.get('rel') in ('stylesheet', 'preload', 'modulepreload'):
                url = attributes.get('href', '')
            if url.startswith(('http://', 'https://', '//')):
                raise SystemExit(f'External rendering dependency in {self.page}')

    class MathMarkup(HTMLParser):
        """Catch TeX accidentally consumed by Markdown, outside literal code."""
        def __init__(self):
            super().__init__()
            self.stack = []

        def handle_starttag(self, tag, attrs):
            if tag not in ('br', 'hr', 'img', 'input', 'meta', 'link', 'wbr', 'source'):
                self.stack.append((tag, dict(attrs)))

        def handle_endtag(self, tag):
            for i in range(len(self.stack) - 1, -1, -1):
                if self.stack[i][0] == tag:
                    del self.stack[i:]
                    break

        def handle_data(self, data):
            if not any(tag == 'article' for tag, _ in self.stack):
                return
            if any(tag in ('code', 'pre', 'script', 'style') or
                   'arithmatex' in attrs.get('class', '').split()
                   for tag, attrs in self.stack):
                return
            if re.search(r'\\(?:frac|sum|prod|theta|epsilon|varepsilon|int|log|operatorname)\b', data):
                raise SystemExit(f'Unparsed LaTeX outside math/code in {self.page}: {data[:100]!r}')

    for path in SITE.rglob('*.html'):
        parser = Assets()
        parser.page = path.relative_to(SITE)
        parser.feed(path.read_text())
        math = MathMarkup()
        math.page = path.relative_to(SITE)
        math.feed(path.read_text())
    print('Generated references and offline rendering asset checks passed.')


if __name__=='__main__':
    main()
