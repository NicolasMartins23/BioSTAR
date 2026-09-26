<template>
  <q-layout view="lHh Lpr lFf">
    <q-header bordered class="bg-dark text-white">
      <q-toolbar>
        <q-btn flat dense round icon="menu" @click="drawerOpen = !drawerOpen" />
        <q-toolbar-title class="text-weight-bold">BioSTAR</q-toolbar-title>
        <q-btn flat round :icon="isDark ? 'light_mode' : 'dark_mode'" @click="toggleDarkMode" />
      </q-toolbar>
    </q-header>
    <q-drawer v-model="drawerOpen" show-if-above bordered :width="240">
      <q-list padding>
        <q-item-label header class="text-uppercase text-weight-bold">Analysis</q-item-label>
        <q-item clickable v-ripple to="/" exact><q-item-section avatar><q-icon name="dashboard" /></q-item-section><q-item-section>Dashboard</q-item-section></q-item>
        <q-item clickable v-ripple to="/sequences"><q-item-section avatar><q-icon name="biotech" /></q-item-section><q-item-section>Sequences</q-item-section></q-item>
        <q-item clickable v-ripple to="/proteins"><q-item-section avatar><q-icon name="science" /></q-item-section><q-item-section>Proteins</q-item-section></q-item>
        <q-item clickable v-ripple to="/mutations"><q-item-section avatar><q-icon name="compare_arrows" /></q-item-section><q-item-section>Mutations</q-item-section></q-item>
        <q-separator class="q-my-md" />
        <q-item clickable v-ripple to="/settings"><q-item-section avatar><q-icon name="settings" /></q-item-section><q-item-section>Settings</q-item-section></q-item>
      </q-list>
    </q-drawer>
    <q-page-container><router-view /></q-page-container>
  </q-layout>
</template>

<script setup lang="ts">
import { computed, ref } from "vue";
import { Dark } from "quasar";

const drawerOpen = ref(true);
const isDark = computed<boolean>(() => Dark.isActive);

const toggleDarkMode = (): void => {
  Dark.toggle();
};
</script>
