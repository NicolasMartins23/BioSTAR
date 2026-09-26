<template>
  <q-page padding>
    <div class="page-shell">
      <div class="text-overline text-primary">Bioinformatics Platform</div>
      <div class="text-h3 text-weight-bold q-mb-sm">Welcome to BioSTAR</div>
      <div class="text-body1 text-grey-6 q-mb-md">
        Analyze biological sequences, proteins, and mutations in one workspace.
      </div>

      <q-chip
        icon="cloud"
        :color="apiStatusColor"
        text-color="white"
      >
        API: {{ apiStatusLabel }}
      </q-chip>

      <div class="row q-col-gutter-lg q-mt-md">
        <div v-for="card in cards" :key="card.title" class="col-12 col-md-4">
          <q-card flat bordered class="analysis-card">
            <q-card-section>
              <q-icon :name="card.icon" size="32px" color="primary" />
              <div class="text-h6 q-mt-md">{{ card.title }}</div>
              <div class="text-body2 text-grey-6 q-mt-sm">{{ card.description }}</div>
            </q-card-section>
            <q-card-actions>
              <q-btn flat color="primary" :label="card.action" :to="card.to" />
            </q-card-actions>
          </q-card>
        </div>
      </div>
    </div>
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
  if (apiStatus.value === "connected") {
    return "Connected";
  }

  if (apiStatus.value === "unavailable") {
    return "Unavailable";
  }

  return "Checking...";
});

const apiStatusColor = computed<string>(() => {
  if (apiStatus.value === "connected") {
    return "positive";
  }

  if (apiStatus.value === "unavailable") {
    return "negative";
  }

  return "warning";
});

const cards: AnalysisCard[] = [
  {
    title: "Sequence Analysis",
    description: "Convert and inspect nucleotide and amino acid sequences.",
    icon: "biotech",
    action: "Open sequences",
    to: "/sequences",
  },
  {
    title: "Protein Analysis",
    description: "Explore protein properties and biological information.",
    icon: "science",
    action: "Open proteins",
    to: "/proteins",
  },
  {
    title: "Mutation Analysis",
    description: "Compare sequences and inspect potential mutations.",
    icon: "compare_arrows",
    action: "Open mutations",
    to: "/mutations",
  },
];

onMounted(async (): Promise<void> => {
  try {
    const health = await getHealth();

    if (health.data?.status !== "ok") {
      apiStatus.value = "unavailable";
      return;
    }

    apiStatus.value = "connected";
  } catch {
    apiStatus.value = "unavailable";
  }
});
</script>

<style scoped>
.page-shell {
  max-width: 1200px;
  margin: 0 auto;
}

.analysis-card {
  height: 100%;
}
</style>
