document.addEventListener("DOMContentLoaded", function() {
    var spoilerElements = document.getElementsByClassName("spoiler");
    for (var i = 0; i < spoilerElements.length; i++) {
      spoilerElements[i].addEventListener("click", function() {
        this.classList.toggle("revealed");
      });
    }
  });