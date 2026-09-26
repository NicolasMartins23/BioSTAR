import { createApp } from "vue";
import { Quasar, Dark } from "quasar";
import "quasar/src/css/index.sass";
import "@quasar/extras/material-icons/material-icons.css";

import App from "./App.vue";
import router from "./router";

Dark.set(true);

createApp(App)
  .use(Quasar)
  .use(router)
  .mount("#q-app");
