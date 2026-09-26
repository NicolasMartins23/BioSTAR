import { defineConfig } from "#q-app";

export default defineConfig(() => ({
  css: [],
  extras: ["material-icons"],
  framework: {},
  build: {
    vueRouterMode: "hash",
  },
  devServer: {
    port: 3000,
    open: false,
  },
}));
