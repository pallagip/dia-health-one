(() => {
  const root = document.querySelector("[data-proxy-viewer]");
  if (!root) return;

  const status = root.querySelector("[data-proxy-status]");
  const resetButton = root.querySelector("[data-proxy-reset]");
  const host = root.querySelector("[data-proxy-canvas]");
  const toggleButtons = Array.from(root.querySelectorAll("[data-proxy-toggle]"));

  const assets = {
    macrocycle: { path: "proxy/inputs/macrocycle_fig8i_tricarboxylate.sdf", format: "sdf" },
    glucose: { path: "proxy/inputs/glucose_model_pose.sdf", format: "sdf" },
  };
  const palette = {
    macrocycle: "#6c53a2",
    glucose: "#d5af5b",
    glucoseHeavy: "#c9922f",
  };
  // 3Dmol 2.5.5 models have no setVisible(); an empty style object hides a model.
  const styles = {
    macrocycle: { stick: { color: palette.macrocycle, radius: 0.17 }, sphere: { scale: 0.22, color: palette.macrocycle } },
    glucose: { stick: { color: palette.glucose, radius: 0.2 }, sphere: { scale: 0.28, color: palette.glucoseHeavy } },
  };

  if (!window.$3Dmol || typeof window.$3Dmol.createViewer !== "function") {
    if (status) status.textContent = "3D module unavailable";
    return;
  }

  const viewer = window.$3Dmol.createViewer(host, {
    backgroundColor: "#fbfafc",
    antialias: true,
    disableFog: false,
  });

  const loaded = {};
  const defaultVisible = { macrocycle: true, glucose: true };
  const visible = { ...defaultVisible };

  async function loadAsset(key) {
    if (loaded[key]) return loaded[key];
    const asset = assets[key];
    const response = await fetch(asset.path, { cache: "force-cache" });
    if (!response.ok) throw new Error(`Could not load ${asset.path}`);
    const text = await response.text();
    const model = viewer.addModel(text, asset.format);
    model.setStyle({}, styles[key]);
    loaded[key] = model;
    return model;
  }

  function applyVisibility() {
    Object.entries(loaded).forEach(([key, model]) => {
      model.setStyle({}, visible[key] ? styles[key] : {});
    });
    toggleButtons.forEach((button) => {
      const key = button.dataset.proxyToggle;
      button.setAttribute("aria-pressed", String(Boolean(visible[key])));
    });
    viewer.render();
  }

  async function boot() {
    try {
      await Promise.all(["macrocycle", "glucose"].map(loadAsset));
      viewer.zoomTo();
      viewer.render();
      if (status) {
        status.textContent = "2 proxy inputs loaded · drag to rotate";
        status.dataset.state = "ready";
      }
    } catch (error) {
      if (status) {
        status.textContent = "3D view unavailable · use the download links";
        status.dataset.state = "partial";
      }
    }
  }

  toggleButtons.forEach((button) => {
    button.addEventListener("click", () => {
      const key = button.dataset.proxyToggle;
      visible[key] = !visible[key];
      applyVisibility();
    });
  });

  if (resetButton) {
    resetButton.addEventListener("click", () => {
      Object.assign(visible, defaultVisible);
      viewer.zoomTo();
      applyVisibility();
    });
  }

  boot();
})();
