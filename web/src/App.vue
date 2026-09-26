<template>
  <q-layout view="lHh Lpr lFf">
    <q-drawer v-model="drawerOpen" show-if-above bordered :width="258">
      <q-scroll-area class="fit">
        <div class="q-pa-md">
          <div class="row items-center q-gutter-sm q-pb-lg">
            <q-avatar color="primary" text-color="white" rounded size="40px"><q-icon name="biotech" /></q-avatar>
            <div><div class="brand-title">BioSTAR</div><div class="text-caption text-grey-6 bio-mono">BIOINFORMATICS</div></div>
          </div>
          <q-list padding>
            <q-item clickable :active="tool === 'home'" @click="tool = 'home'" active-class="text-primary bg-grey-2">
              <q-item-section avatar><q-icon name="home" /></q-item-section><q-item-section>Start</q-item-section>
            </q-item>
            <q-item-label header>GENE</q-item-label>
            <q-item v-for="item in geneTools" :key="item.id" clickable :active="tool === item.id" @click="tool = item.id">
              <q-item-section avatar><q-icon :name="item.icon" /></q-item-section><q-item-section>{{ item.label }}</q-item-section>
            </q-item>
            <q-item-label header>PROTEIN</q-item-label>
            <q-item clickable :active="tool === 'protein'" @click="tool = 'protein'">
              <q-item-section avatar><q-icon name="biotech" /></q-item-section><q-item-section>Protein Analysis</q-item-section>
            </q-item>
            <q-item-label header>OTHER</q-item-label>
            <q-item clickable tag="a" href="/docs" target="_blank">
              <q-item-section avatar><q-icon name="menu_book" /></q-item-section><q-item-section>API Documentation</q-item-section>
            </q-item>
          </q-list>
          <q-btn flat no-caps icon="code" label="GitHub" class="full-width q-mt-md" href="https://github.com/NicolasMartins23/BioSTAR" target="_blank" />
        </div>
      </q-scroll-area>
    </q-drawer>

    <q-header bordered class="bg-background text-grey-7">
      <q-toolbar>
        <q-btn flat round dense icon="menu" @click="drawerOpen = !drawerOpen" />
        <q-toolbar-title class="bio-mono text-caption">BioSTAR / {{ toolTitle }}</q-toolbar-title>
        <q-btn flat round dense :icon="$q.dark.isActive ? 'light_mode' : 'dark_mode'" @click="$q.dark.toggle()" />
      </q-toolbar>
    </q-header>

    <q-page-container><q-page class="q-pa-md q-pa-lg-xl"><div class="page-width">
      <HomePage v-if="tool === 'home'" @select="tool = $event" />
      <ConversionPage v-else-if="isConversion" :tool="tool" />
      <ProteinPage v-else-if="tool === 'protein'" />
      <MutationPage v-else />
    </div></q-page></q-page-container>
  </q-layout>
</template>

<script setup lang="ts">
import { computed, ref } from "vue";
import HomePage from "./components/HomePage.vue";
import ConversionPage from "./components/ConversionPage.vue";
import ProteinPage from "./components/ProteinPage.vue";
import MutationPage from "./components/MutationPage.vue";

type Tool = "home" | "dna-rna" | "dna-protein" | "rna-protein" | "rna-dna" | "protein" | "mutation";
const drawerOpen = ref(true);
const tool = ref<Tool>("home");
const geneTools: Array<{ id: Exclude<Tool, "home"|"protein"|"mutation"> | "mutation"; label: string; icon: string }> = [
  { id: "dna-protein", label: "DNA → Protein", icon: "science" },
  { id: "dna-rna", label: "DNA → RNA", icon: "swap_horiz" },
  { id: "rna-protein", label: "RNA → Protein", icon: "swap_horiz" },
  { id: "rna-dna", label: "RNA → DNA", icon: "swap_horiz" },
  { id: "mutation", label: "Mutation Compare", icon: "compare_arrows" },
];
const isConversion = computed(() => ["dna-rna","dna-protein","rna-protein","rna-dna"].includes(tool.value));
const toolTitle = computed(() => ({home:"start","dna-rna":"dna → rna","dna-protein":"dna → protein","rna-protein":"rna → protein","rna-dna":"rna → dna",protein:"protein analysis",mutation:"mutation compare"}[tool.value]));
</script>

<style scoped>.page-width{width:min(1250px,100%);margin:0 auto}</style>
