// "Your values" form of guides/rclone_ovh_s3.md.
// The reader picks a system (Linux, Mac, Windows); blocks marked data-only="..." are
// shown for that system only, and code blocks marked .rclone-tpl and [data-rclone]
// spans are rewritten with the reader's values. Values are kept in localStorage,
// never sent anywhere. Initialised through Material's document$ observable, so it
// also works with navigation.instant (the page body is replaced without a reload).
(function () {
  "use strict";

  const KEY = "startbio-rclone-ovh-v2";
  const OSES = ["linux", "mac", "win"];
  // os: systems that show the field; d: default (an object when it differs by system).
  const FIELDS = [
    { k: "remote", l: "rclone configuration name", os: OSES, d: "ovh", h: "Short name used in all commands (step 2)." },
    { k: "bucket", l: "Bucket name", os: OSES, d: "my-bucket", ex: true, h: "Provided with your keys." },
    { k: "region", l: "Bucket region", os: OSES, d: "eu-west-par", h: "Paris 3-AZ: the only region with Cold Archive." },
    { k: "src", l: "Source directory (real data)", os: OSES, d: { linux: "/path/to/data", mac: "/path/to/data", win: "D:\\path\\to\\data" }, ex: true, h: "Full path on this machine (step 6)." },
    { k: "dst", l: "Destination directory in the bucket (real data)", os: OSES, d: "my-dataset", ex: true, h: "Without a leading / (step 6)." },
    { k: "day", l: "Bandwidth on weekdays, 8 am – 8 pm", os: OSES, d: "50M", h: "In bytes per second: 50M = 50 MiB/s." },
    { k: "night", l: "Bandwidth at night and at weekends", os: OSES, d: "80M", h: "Lower it if the network suffers." }
  ];
  const RATE = /^\d+(\.\d+)?[KMG]?$/;
  const PLACEHOLDER = /\{\{(\w+)\}\}/g;

  function load() {
    let s = {};
    try { s = JSON.parse(localStorage.getItem(KEY) || "{}") || {}; } catch (e) { s = {}; }
    if (!OSES.includes(s.os)) delete s.os;
    return s;
  }
  function save(s) { try { localStorage.setItem(KEY, JSON.stringify(s)); } catch (e) { /* storage unavailable */ } }

  const def = (f, os) => (typeof f.d === "string" ? f.d : f.d[os]);
  // Fields whose default depends on the system are stored per system.
  const skey = (f, os) => (typeof f.d === "string" ? f.k : os + "." + f.k);
  const raw = (s, f, os) => (typeof s[skey(f, os)] === "string" ? s[skey(f, os)] : def(f, os));

  // ~ becomes $HOME (Windows: $env:USERPROFILE), which still expands inside double
  // quotes; no trailing separator.
  function path(os, s) {
    s = s.trim();
    if (os === "win") return s.replace(/^~(?=[\\/]|$)/, "$env:USERPROFILE").replace(/[\\/]+$/, "");
    return s.replace(/^~(?=\/|$)/, "$HOME").replace(/\/+$/, "");
  }

  function values(s, os) {
    const v = {};
    FIELDS.forEach(f => { if (f.os.includes(os)) v[f.k] = raw(s, f, os).trim(); });
    v.src = path(os, v.src);
    v.dst = v.dst.replace(/^[\\/]+|[\\/]+$/g, "");
    // Exercise data and logs always go to the home folder: no field for them.
    const sep = os === "win" ? "\\" : "/";
    const home = os === "win" ? "$env:USERPROFILE" : "$HOME";
    v.exo = home + sep + "rclone-exercise";
    v.logs = home + sep + "rclone-logs";
    v.logcopy = v.logs + sep + "exercise_copy.log";
    v.newfile = home + sep + "new_sample.bam"; // step 5
    v.newbatch = home + sep + "new_batch";
    v.day = RATE.test(v.day) ? v.day : "50M";
    v.night = RATE.test(v.night) ? v.night : "80M";
    v.bw = ["Mon", "Tue", "Wed", "Thu", "Fri"]
      .map(d => d + "-08:00," + v.day + " " + d + "-20:00," + v.night).join(" ");
    v.endpointhost = "s3." + (v.region || "eu-west-par") + ".io.cloud.ovh.net";
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
    // The attribute that drives which blocks are shown (see rclone-form.css).
    const scope = form.closest("article") || document.body;

    const s = load();
    const codes = Array.from(document.querySelectorAll(".rclone-tpl code"));
    codes.forEach(c => { c.dataset.tpl = c.textContent; });

    function build() {
      box.textContent = "";
      if (!s.os) return;
      const os = s.os;
      FIELDS.filter(f => f.os.includes(os)).forEach(f => {
        const w = document.createElement("div");
        w.className = "rclone-field";
        const lab = document.createElement("label");
        const inp = document.createElement("input");
        const hint = document.createElement("div");
        inp.type = "text"; inp.spellcheck = false; inp.autocomplete = "off";
        inp.id = "rclone-f-" + f.k; lab.htmlFor = inp.id;
        lab.textContent = f.l; inp.value = raw(s, f, os);
        hint.className = "rclone-hint";
        w.append(lab, inp, hint);
        box.appendChild(w);
        f.inp = inp; f.hint = hint;
        inp.addEventListener("input", () => { s[skey(f, os)] = inp.value; save(s); update(); });
      });
    }

    function update() {
      form.querySelectorAll("[data-os-pick]").forEach(b => b.setAttribute("aria-pressed", String(b.dataset.osPick === s.os)));
      scope.dataset.rcloneOs = s.os || "none";
      // Before a system is chosen, only the blocks common to all systems are visible;
      // they use shared fields only, so the Linux defaults fill them correctly.
      const os = s.os || "linux", v = values(s, os);
      if (s.os) FIELDS.filter(f => f.os.includes(os) && f.inp).forEach(f => {
        const r = raw(s, f, os).trim();
        let bad = false, msg = f.h;
        if ((f.k === "day" || f.k === "night") && !RATE.test(r)) { bad = true; msg = "Expected format: 50M, 80M, 500K…"; }
        if (!r) { bad = true; msg = "This field must not be empty."; }
        const ex = !!f.ex && r === def(f, os);
        f.inp.classList.toggle("rclone-bad", bad);
        f.inp.classList.toggle("rclone-ex", ex && !bad);
        f.hint.textContent = ex && !bad ? "Example value, to replace. " + msg : msg;
        f.hint.classList.toggle("rclone-warn", bad || ex);
      });
      codes.forEach(c => render(c, c.dataset.tpl, v));
      document.querySelectorAll("[data-rclone]").forEach(e => { e.textContent = v[e.dataset.rclone] ?? ""; });
    }

    form.querySelectorAll("[data-os-pick]").forEach(b => b.addEventListener("click", () => {
      s.os = b.dataset.osPick; save(s); build(); update();
    }));
    form.querySelector(".rclone-reset").addEventListener("click", () => {
      if (!s.os) return;
      FIELDS.filter(f => f.os.includes(s.os)).forEach(f => { delete s[skey(f, s.os)]; });
      save(s); build(); update();
    });
    build();
    update();
  }

  if (typeof document$ !== "undefined") document$.subscribe(init);
  else document.addEventListener("DOMContentLoaded", init);
})();
