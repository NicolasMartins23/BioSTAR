import { createApp } from "vue";
import { Quasar } from "quasar";
import "quasar/src/css/index.sass";
import "@quasar/extras/material-icons/material-icons.css";

import App from "./App.vue";

createApp(App)
  .use(Quasar)
  .mount("#q-app");
