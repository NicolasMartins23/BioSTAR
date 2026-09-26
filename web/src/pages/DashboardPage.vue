<template>
  <q-page class="biostar-page">
    <main class="biostar-content">
      <section>
        <div class="biostar-eyebrow">Bioinformatics Analysis Suite</div>
        <h1 class="biostar-title">Biological analysis, in one workspace.</h1>
        <p class="biostar-subtitle">
          BioSTAR provides sequence conversion and biochemical analysis tools
          through a unified scientific workspace.
        </p>
      </section>

      <section class="q-mt-xl">
        <div class="biostar-panel">
          <div class="biostar-panel__body row items-center justify-between q-col-gutter-lg">
            <div class="col-12 col-sm">
              <div class="biostar-section-label">System status</div>
              <div class="text-body2 text-grey-6 q-mt-xs">
                Connection to the BioSTAR analysis API
              </div>
            </div>
            <div class="col-auto">
              <q-chip square :color="apiStatusColor" text-color="white" :icon="apiStatusIcon">
                {{ apiStatusLabel }}
              </q-chip>
            </div>
          </div>
        </div>
      </section>

      <section class="q-mt-xl">
        <div class="biostar-section-label q-mb-md">Analysis tools</div>
        <div class="row q-col-gutter-md">
          <div v-for="card in cards" :key="card.title" class="col-12 col-md-4">
            <q-card flat class="biostar-panel tool-card">
              <q-card-section>
                <q-icon :name="card.icon" size="30px" color="primary" />
                <div class="text-h6 text-weight-bold q-mt-lg">{{ card.title }}</div>
                <div class="text-body2 text-grey-6 q-mt-sm tool-card__description">
                  {{ card.description }}
                </div>
              </q-card-section>
              <q-card-actions class="q-px-md q-pb-md">
                <q-btn flat no-caps color="primary" :label="card.action" :to="card.to" />
              </q-card-actions>
            </q-card>
          </div>
        </div>
      </section>
    </main>
  </q-page>
</template>

<script setup lang="ts">
import { computed, onMounted, ref } from "vue";
import { getHealth } from "../services/api";

type ApiStatus = "checking" | "connected" | "unavailable";

interface AnalysisCard {
  title: string;
  description: string;
  icon: string;
  action: string;
  to: string;
}

const apiStatus = ref<ApiStatus>("checking");

const apiStatusLabel = computed<string>(() => {
  if (apiStatus.value === "connected") return "API connected";
  if (apiStatus.value === "unavailable") return "API unavailable";
  return "Checking API";
});

const apiStatusColor = computed<string>(() => {
  if (apiStatus.value === "connected") return "positive";
  if (apiStatus.value === "unavailable") return "negative";
  return "warning";
});

const apiStatusIcon = computed<string>(() => {
  if (apiStatus.value === "connected") return "check_circle";
  if (apiStatus.value === "unavailable") return "error";
  return "sync";
});

const cards: AnalysisCard[] = [
  {
    title: "Sequence conversion",
    description: "Convert DNA and RNA sequences and translate nucleotide sequences into proteins.",
    icon: "biotech",
    action: "Open analysis",
    to: "/sequences",
  },
  {
    title: "Protein analysis",
    description: "Calculate biochemical properties and inspect protein composition.",
    icon: "science",
    action: "Open analysis",
    to: "/proteins",
  },
  {
    title: "Mutation analysis",
    description: "Compare biological sequences and identify sequence-level changes.",
    icon: "compare_arrows",
    action: "Open analysis",
    to: "/mutations",
  },
];

onMounted(async (): Promise<void> => {
  try {
    const health = await getHealth();
    apiStatus.value = health.data?.status === "ok" ? "connected" : "unavailable";
  } catch {
    apiStatus.value = "unavailable";
  }
});
</script>

<style scoped>
.tool-card {
  height: 100%;
  transition: border-color 120ms ease, transform 120ms ease;
}

.tool-card:hover {
  border-color: rgb(23 107 135 / 35%);
  transform: translateY(-1px);
}

.tool-card__description {
  min-height: 48px;
  line-height: 1.55;
}
</style>
