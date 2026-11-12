(function () {
  function renderHypraMath() {
    if (typeof katex === 'undefined') return;
    document.querySelectorAll('.arithmatex').forEach(function (element) {
      if (element.dataset.hypraMathRendered) return;
      var text = element.textContent.trim();
      var display = text.startsWith('\\[') && text.endsWith('\\]');
      var inline = text.startsWith('\\(') && text.endsWith('\\)');
      if (!display && !inline) return;
      katex.render(text.slice(2, -2), element, {
        displayMode: display,
        throwOnError: false,
        strict: 'warn'
      });
      element.dataset.hypraMathRendered = 'true';
    });
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', renderHypraMath);
  } else {
    renderHypraMath();
  }
  if (typeof document$ !== 'undefined') document$.subscribe(renderHypraMath);
})();
