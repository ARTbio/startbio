// "Your values" form of guides/rclone_ovh_s3.md.
// Rewrites the code blocks marked .rclone-tpl and the [data-rclone] spans with the
// values typed by the reader. Values are kept in localStorage, never sent anywhere.
// Initialised through Material's document$ observable, so it also works with
// navigation.instant (the page body is replaced without a full reload).
(function () {
  "use strict";

  const KEY = "startbio-rclone-ovh-v1";
  const FIELDS = [
    { k: "remote", l: "rclone remote name", d: "ovh", h: "Short name given to the configuration (step 2)." },
    { k: "bucket", l: "Bucket name", d: "my-bucket", ex: true, h: "Provided with your keys." },
    { k: "region", l: "Bucket region", d: "eu-west-par", h: "Paris 3-AZ: the only region with Cold Archive." },
    { k: "disk", l: "Data disk for the exercise", d: "/data_disk", ex: true, h: "A directory you can write to, on a data disk." },
    { k: "vol", l: "Exercise volume (GB)", d: "100", h: "100 by default; 500 for a more thorough exercise." },
    { k: "src", l: "Directory of the real data", d: "/path/to/data", ex: true, h: "Full path on the server." },
    { k: "dst", l: "Destination directory in the bucket", d: "my-dataset", ex: true, h: "Without a leading /." },
    { k: "day", l: "Bandwidth on weekdays, 8 am – 8 pm", d: "50M", h: "In bytes per second: 50M = 50 MiB/s." },
    { k: "night", l: "Bandwidth at night and at weekends", d: "80M", h: "Lower it if the site network suffers." },
    { k: "logs", l: "Log directory", d: "$HOME/rclone-logs", h: "Created at step 4 b)." }
  ];
  const RATE = /^\d+(\.\d+)?[KMG]?$/;
  const PLACEHOLDER = /\{\{(\w+)\}\}/g;

  function load() {
    let s = {};
    try { s = JSON.parse(localStorage.getItem(KEY) || "{}") || {}; } catch (e) { s = {}; }
    FIELDS.forEach(f => { if (typeof s[f.k] !== "string") s[f.k] = f.d; });
    return s;
  }
  function save(s) { try { localStorage.setItem(KEY, JSON.stringify(s)); } catch (e) { /* storage unavailable */ } }

  // ~ becomes $HOME (which still expands inside double quotes); no trailing slash.
  const home = s => s.trim().replace(/^~(?=\/|$)/, "$HOME").replace(/\/+$/, "");
  const NORM = {
    logs: s => home(s) || "$HOME/rclone-logs",
    disk: home,
    src: home,
    dst: s => s.trim().replace(/^\/+|\/+$/g, "")
  };

  function values(s) {
    const v = {};
    FIELDS.forEach(f => { v[f.k] = NORM[f.k] ? NORM[f.k](s[f.k]) : s[f.k].trim(); });
    const n = parseInt(v.vol, 10);
    v.vol = n > 0 ? String(n) : "100";
    v.need = String(Math.ceil((+v.vol < 10 ? 2 : +v.vol) * 11 / 10));
    v.exo = (v.disk || "/data_disk") + "/rclone-exercise";
    v.day = RATE.test(v.day) ? v.day : "50M";
    v.night = RATE.test(v.night) ? v.night : "80M";
    v.bw = ["Mon", "Tue", "Wed", "Thu", "Fri"]
      .map(d => d + "-08:00," + v.day + " " + d + "-20:00," + v.night).join(" ");
    return v;
  }

  // Rebuilds a code element from its template, one span per line (Material lays out
  // code lines as grid rows); values are text nodes in an inline span.
  function render(code, tpl, v) {
    code.textContent = "";
    tpl.split(/(?<=\n)/).forEach(line => {
      const row = document.createElement("span");
      let last = 0, m;
      PLACEHOLDER.lastIndex = 0;
      while ((m = PLACEHOLDER.exec(line))) {
        if (m.index > last) row.appendChild(document.createTextNode(line.slice(last, m.index)));
        const val = document.createElement("span");
        val.className = "rclone-v";
        val.textContent = v[m[1]] ?? "";
        row.appendChild(val);
        last = PLACEHOLDER.lastIndex;
      }
      row.appendChild(document.createTextNode(line.slice(last)));
      code.appendChild(row);
    });
  }

  function init() {
    const form = document.getElementById("rclone-form");
    if (!form) return;
    const box = form.querySelector(".rclone-fields");
    if (!box || box.dataset.ready) return;
    box.dataset.ready = "1";

    const s = load();
    const codes = Array.from(document.querySelectorAll(".rclone-tpl code"));
    codes.forEach(c => { c.dataset.tpl = c.textContent; });

    FIELDS.forEach(f => {
      const w = document.createElement("div");
      w.className = "rclone-field";
      const lab = document.createElement("label");
      const inp = document.createElement("input");
      const hint = document.createElement("div");
      inp.type = "text"; inp.spellcheck = false; inp.autocomplete = "off";
      inp.id = "rclone-f-" + f.k; lab.htmlFor = inp.id;
      lab.textContent = f.l; inp.value = s[f.k];
      hint.className = "rclone-hint";
      w.append(lab, inp, hint);
      box.appendChild(w);
      f.inp = inp; f.hint = hint;
      inp.addEventListener("input", () => { s[f.k] = inp.value; save(s); update(); });
    });

    function update() {
      const v = values(s);
      FIELDS.forEach(f => {
        const raw = s[f.k].trim();
        let bad = false, msg = f.h;
        if ((f.k === "day" || f.k === "night") && !RATE.test(raw)) { bad = true; msg = "Expected format: 50M, 80M, 500K…"; }
        if (f.k === "vol" && !(parseInt(raw, 10) > 0)) { bad = true; msg = "A whole number of GB, for instance 100."; }
        if (!raw && f.k !== "dst") { bad = true; msg = "This field must not be empty."; }
        const homeDisk = f.k === "disk" && v.disk.startsWith("$HOME");
        if (homeDisk) msg = "Your home directory is often on the system disk: avoid it.";
        const ex = !!f.ex && raw === f.d;
        f.inp.classList.toggle("rclone-bad", bad);
        f.inp.classList.toggle("rclone-ex", ex && !bad);
        f.hint.textContent = ex && !bad ? "Example value, to replace. " + msg : msg;
        f.hint.classList.toggle("rclone-warn", bad || ex || homeDisk);
      });
      codes.forEach(c => render(c, c.dataset.tpl, v));
      document.querySelectorAll("[data-rclone]").forEach(e => { e.textContent = v[e.dataset.rclone] ?? ""; });
    }

    form.querySelector(".rclone-reset").addEventListener("click", () => {
      FIELDS.forEach(f => { s[f.k] = f.d; f.inp.value = f.d; });
      save(s); update();
    });
    update();
  }

  if (typeof document$ !== "undefined") document$.subscribe(init);
  else document.addEventListener("DOMContentLoaded", init);
})();
