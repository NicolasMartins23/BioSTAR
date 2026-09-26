import { createApp } from "vue";
import { Quasar, Dark } from "quasar";
import "quasar/src/css/index.sass";
import "@quasar/extras/material-icons/material-icons.css";
import "./css/biostar.scss";

import App from "./App.vue";
import router from "./router";

Dark.set(false);

createApp(App)
  .use(Quasar, {
    config: {
      brand: {
        primary: "#176b87",
        secondary: "#238b8f",
        accent: "#4f8a5b",
        dark: "#12324a",
        positive: "#4f8a5b",
        negative: "#b94a48",
        info: "#176b87",
        warning: "#b5792d",
      },
    },
  })
  .use(router)
  .mount("#q-app");
