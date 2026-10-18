(() => {
  const root = document.querySelector("[data-structure-viewer]");
  if (!root) return;

  const status = root.querySelector("[data-viewer-status]");
  const resetButton = root.querySelector("[data-viewer-reset]");
  const palette = {
    insulinA: "#9180ab",
    insulinB: "#b4a7c7",
    macrocycle: "#6c53a2",
    glucoside: "#92ae6e",
    glucose: "#d5af5b",
  };
  const panels = [
    { id: "nnc-open-viewer", path: "structures/nnc2215-open-model.pdb" },
    { id: "nnc-closed-viewer", path: "structures/nnc2215-closed-model.pdb" },
  ];
  const viewers = [];

  if (!window.$3Dmol || typeof window.$3Dmol.createViewer !== "function") {
    status.textContent = "3D module unavailable";
    return;
  }

  async function loadPanel(panel) {
    const host = document.getElementById(panel.id);
    const response = await fetch(panel.path, { cache: "force-cache" });
    if (!response.ok) throw new Error(`Could not load ${panel.path}`);
    const pdb = await response.text();
    const viewer = window.$3Dmol.createViewer(host, {
      backgroundColor: "#fbfafc",
      antialias: true,
      disableFog: false,
    });

    viewer.addModel(pdb, "pdb");
    viewer.setStyle({ chain: "A" }, { cartoon: { color: palette.insulinA } });
    viewer.setStyle({ chain: "B" }, { cartoon: { color: palette.insulinB } });
    viewer.setStyle({ chain: "B", resn: "PHX" }, { stick: { color: palette.glucoside, radius: 0.17 } });
    viewer.setStyle({ chain: "B", resn: "LYX" }, { stick: { color: palette.macrocycle, radius: 0.15 } });
    viewer.setStyle({ chain: "C" }, { stick: { color: palette.macrocycle, radius: 0.16 } });
    viewer.setStyle({ chain: "D" }, { stick: { color: palette.glucose, radius: 0.18 } });
    viewer.zoomTo();
    viewer.render();
    viewers.push(viewer);
    return viewer;
  }

  Promise.allSettled(panels.map(loadPanel)).then((results) => {
    const ready = results.filter((result) => result.status === "fulfilled").length;
    status.textContent = ready === panels.length
      ? "2 source models loaded · rotate to inspect"
      : ready === 1
        ? "1 model loaded · use the source links if needed"
        : "3D view unavailable · use the source links";
    status.dataset.state = ready === panels.length ? "ready" : "partial";
  });

  resetButton.addEventListener("click", () => {
    viewers.forEach((viewer) => {
      viewer.zoomTo();
      viewer.render();
    });
  });
})();
