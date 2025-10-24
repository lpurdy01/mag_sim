(function () {
  const scheduleTypeset = () => {
    const tryTypeset = () => {
      if (window.MathJax && typeof window.MathJax.typesetPromise === 'function') {
        window.MathJax.typesetPromise().catch((err) => console.error('MathJax typeset failed', err));
      } else {
        window.setTimeout(tryTypeset, 50);
      }
    };
    tryTypeset();
  };

  if (window.document$ && typeof window.document$.subscribe === 'function') {
    window.document$.subscribe(scheduleTypeset);
  } else {
    document.addEventListener('DOMContentLoaded', scheduleTypeset);
  }

  window.addEventListener('load', scheduleTypeset);
})();
