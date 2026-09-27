<template>
  <q-page class="biostar-page">
    <div class="biostar-content">
      <h1 class="biostar-title">BioSTAR</h1>
      <p class="biostar-subtitle">Scientific workspace</p>
    </div>
  </q-page>
</template>

<script setup lang="ts">
import { onMounted, ref } from "vue";
import { getHealth } from "../services/api";

const apiConnected = ref<boolean>(false);

onMounted(async (): Promise<void> => {
  try {
    const response = await getHealth();
    apiConnected.value = response.data?.status === "ok";
  } catch {
    apiConnected.value = false;
  }
});
</script>
