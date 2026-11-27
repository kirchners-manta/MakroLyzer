/* Enhance ordinary Sphinx figures; all drawings remain readable without JS. */
(() => {
  "use strict";

  function initializeGalleries() {
    document.querySelectorAll(".linker-gallery").forEach((gallery) => {
      const slides = Array.from(gallery.querySelectorAll(".linker-slide"));
      const controls = gallery.querySelector(".linker-gallery-controls");
      if (!controls || slides.length < 2) return;

      const status = controls.querySelector(".linker-current");
      let current = 0;
      gallery.setAttribute("role", "region");
      gallery.setAttribute("aria-roledescription", "carousel");
      gallery.setAttribute("aria-label", "Lignin linker gallery");
      gallery.tabIndex = 0;

      // Wrap at either end so readers can browse the complete motif set.
      function show(index) {
        current = (index + slides.length) % slides.length;
        slides.forEach((slide, position) => {
          slide.hidden = position !== current;
          slide.setAttribute("role", "group");
          slide.setAttribute("aria-roledescription", "slide");
          slide.setAttribute("aria-label", `${position + 1} of ${slides.length}`);
        });
        const caption = slides[current].querySelector(".caption-text");
        status.textContent = caption.textContent.trim();
      }

      controls.querySelector(".linker-previous").addEventListener("click", () => show(current - 1));
      controls.querySelector(".linker-next").addEventListener("click", () => show(current + 1));
      gallery.addEventListener("keydown", (event) => {
        const destinations = {
          ArrowLeft: current - 1,
          ArrowRight: current + 1,
          Home: 0,
          End: slides.length - 1,
        };
        if (Object.hasOwn(destinations, event.key)) {
          event.preventDefault();
          show(destinations[event.key]);
        }
      });
      show(0);
      controls.hidden = false;
      gallery.classList.add("is-ready");
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initializeGalleries);
  } else {
    initializeGalleries();
  }
})();
