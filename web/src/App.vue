<template>
  <q-layout view="hHh Lpr lFf">
    <q-header class="biostar-header">
      <q-toolbar class="biostar-toolbar">
        <q-btn
          flat
          dense
          round
          icon="menu"
          aria-label="Toggle navigation"
          class="q-mr-sm"
          @click="drawerOpen = !drawerOpen"
        />

        <q-toolbar-title class="biostar-brand">
          <div class="biostar-brand__mark">B</div>
          <div>
            <div class="biostar-brand__name">BioSTAR</div>
            <div class="biostar-brand__subtitle">Bioinformatics Analysis Suite</div>
          </div>
        </q-toolbar-title>

        <div class="biostar-header__status gt-xs">
          <span class="status-dot" />
          Analysis environment
        </div>

        <q-btn
          flat
          round
          :icon="isDark ? 'light_mode' : 'dark_mode'"
          :aria-label="isDark ? 'Use light theme' : 'Use dark theme'"
          @click="toggleDarkMode"
        >
          <q-tooltip>{{ isDark ? "Light theme" : "Dark theme" }}</q-tooltip>
        </q-btn>
      </q-toolbar>
    </q-header>

    <q-drawer v-model="drawerOpen" show-if-above bordered :width="250" class="biostar-drawer">
      <div class="biostar-drawer__intro">
        <div class="biostar-section-label">Workspace</div>
        <div class="text-caption text-grey-6 q-mt-xs">
          Biological sequence analysis
        </div>
      </div>

      <q-list padding>
        <q-item-label header class="biostar-nav-label">Analysis</q-item-label>

        <q-item clickable v-ripple to="/" exact class="biostar-nav-item">
          <q-item-section avatar><q-icon name="dashboard" /></q-item-section>
          <q-item-section>Overview</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/sequences" class="biostar-nav-item">
          <q-item-section avatar><q-icon name="biotech" /></q-item-section>
          <q-item-section>Sequence conversion</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/proteins" class="biostar-nav-item">
          <q-item-section avatar><q-icon name="science" /></q-item-section>
          <q-item-section>Protein analysis</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/mutations" class="biostar-nav-item">
          <q-item-section avatar><q-icon name="compare_arrows" /></q-item-section>
          <q-item-section>Mutation analysis</q-item-section>
        </q-item>

        <q-separator class="q-my-md" />

        <q-item-label header class="biostar-nav-label">System</q-item-label>
        <q-item clickable v-ripple to="/settings" class="biostar-nav-item">
          <q-item-section avatar><q-icon name="settings" /></q-item-section>
          <q-item-section>Settings</q-item-section>
        </q-item>
      </q-list>

      <div class="biostar-drawer__footer">
        <div class="text-caption text-grey-6">BioSTAR API</div>
        <div class="text-caption text-weight-medium">Scientific analysis workspace</div>
      </div>
    </q-drawer>

    <q-page-container>
      <router-view />
    </q-page-container>
  </q-layout>
</template>

<script setup lang="ts">
import { computed, ref } from "vue";
import { Dark } from "quasar";

const drawerOpen = ref<boolean>(true);
const isDark = computed<boolean>(() => Dark.isActive);

const toggleDarkMode = (): void => {
  Dark.toggle();
};
</script>

<style scoped>
.biostar-header {
  background: var(--biostar-surface);
  color: var(--biostar-text);
  border-bottom: 1px solid var(--biostar-border);
}

.biostar-toolbar {
  min-height: 64px;
  max-width: 1440px;
  margin: 0 auto;
}

.biostar-brand {
  display: flex;
  align-items: center;
  gap: 11px;
}

.biostar-brand__mark {
  display: grid;
  width: 34px;
  height: 34px;
  place-items: center;
  border-radius: 4px;
  background: var(--q-primary);
  color: white;
  font-size: 18px;
  font-weight: 800;
}

.biostar-brand__name {
  font-size: 1.05rem;
  font-weight: 800;
  letter-spacing: 0.04em;
}

.biostar-brand__subtitle {
  color: var(--biostar-muted);
  font-size: 0.68rem;
  line-height: 1.2;
}

.biostar-header__status {
  margin-right: 20px;
  color: var(--biostar-muted);
  font-size: 0.75rem;
}

.status-dot {
  display: inline-block;
  width: 7px;
  height: 7px;
  margin-right: 6px;
  border-radius: 50%;
  background: var(--biostar-green);
}

.biostar-drawer {
  background: var(--biostar-surface);
}

.biostar-drawer__intro {
  padding: 24px 20px 12px;
}

.biostar-nav-label {
  color: var(--biostar-muted);
  font-size: 0.68rem;
  font-weight: 700;
  letter-spacing: 0.12em;
  text-transform: uppercase;
}

.biostar-nav-item {
  min-height: 42px;
  margin: 2px 8px;
  border-radius: 4px;
  color: var(--biostar-text);
}

.biostar-nav-item.q-router-link--active {
  background: rgb(23 107 135 / 9%);
  color: var(--q-primary);
  font-weight: 600;
}

.biostar-drawer__footer {
  position: absolute;
  right: 20px;
  bottom: 18px;
  left: 20px;
  padding-top: 14px;
  border-top: 1px solid var(--biostar-border);
}
</style>
