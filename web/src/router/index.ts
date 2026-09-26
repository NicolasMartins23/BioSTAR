import { createRouter, createWebHashHistory } from "vue-router";

import DashboardPage from "../pages/DashboardPage.vue";
import MutationsPage from "../pages/MutationsPage.vue";
import ProteinsPage from "../pages/ProteinsPage.vue";
import SequencesPage from "../pages/SequencesPage.vue";
import SettingsPage from "../pages/SettingsPage.vue";

const router = createRouter({
  history: createWebHashHistory(),
  routes: [
    { path: "/", component: DashboardPage },
    { path: "/sequences", component: SequencesPage },
    { path: "/proteins", component: ProteinsPage },
    { path: "/mutations", component: MutationsPage },
    { path: "/settings", component: SettingsPage },
  ],
});

export default router;
