import { defineConfig } from "#q-app";

export default defineConfig(() => ({
  boot: [],
  css: ["app.scss"],
  extras: ["material-icons", "mdi-v7"],
  framework: {
    plugins: ["Notify"],
  },
  build: {
    vueRouterMode: "hash",
  },
  devServer: {
    port: 3000,
    open: false,
  },
  htmlVariables: {
    title: "BioSTAR",
  },
}));
