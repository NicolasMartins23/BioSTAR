import { createApp } from "vue";
import { Quasar, Notify } from "quasar";
import "quasar/src/css/index.sass";
import "@quasar/extras/material-icons/material-icons.css";
import App from "./App.vue";
import "./css/app.scss";

createApp(App).use(Quasar, { plugins: { Notify } }).mount("#q-app");
