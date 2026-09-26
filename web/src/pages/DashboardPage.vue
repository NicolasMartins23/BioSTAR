<template>
  <q-page class="biostar-page">
    <main class="biostar-content">
      <section class="hero">
        <div class="hero-copy">
          <div class="biostar-eyebrow">BioSTAR · Bioinformatics resource</div>
          <h1 class="biostar-title">Explore biological data.<br />Run the analysis.</h1>
          <p class="biostar-subtitle">
            A focused workspace for sequence conversion and biochemical protein analysis,
            built around the BioSTAR scientific computing engine.
          </p>
        </div>

        <div class="hero-art" aria-hidden="true">
          <div class="helix">
            <i v-for="n in 7" :key="n" :style="{ '--n': n }"></i>
          </div>
        </div>
      </section>

      <section class="status-line">
        <span class="status-marker" :class="apiStatus"></span>
        <span class="status-label">BioSTAR API</span>
        <span class="status-value">{{ apiStatusLabel }}</span>
      </section>

      <section class="tools">
        <div class="tools-heading">
          <div>
            <div class="biostar-eyebrow">Analysis tools</div>
            <h2>Choose a workflow</h2>
          </div>
        </div>

        <div class="tool-list">
          <q-item
            v-for="(card, index) in cards"
            :key="card.title"
            clickable
            v-ripple
            :to="card.to"
            class="tool-row"
          >
            <q-item-section side class="tool-index">0{{ index + 1 }}</q-item-section>
            <q-item-section avatar>
              <q-icon :name="card.icon" size="28px" color="primary" />
            </q-item-section>
            <q-item-section>
              <q-item-label class="tool-title">{{ card.title }}</q-item-label>
              <q-item-label caption>{{ card.description }}</q-item-label>
            </q-item-section>
            <q-item-section side>
              <q-icon name="arrow_forward" color="primary" />
            </q-item-section>
          </q-item>
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
  to: string;
}

const apiStatus = ref<ApiStatus>("checking");

const apiStatusLabel = computed<string>(() => {
  if (apiStatus.value === "connected") return "Connected";
  if (apiStatus.value === "unavailable") return "Unavailable";
  return "Checking";
});

const cards: AnalysisCard[] = [
  {
    title: "Sequence conversion",
    description: "DNA ↔ RNA conversion and nucleotide-to-protein translation.",
    icon: "biotech",
    to: "/sequences",
  },
  {
    title: "Protein analysis",
    description: "Physicochemical properties, composition and biochemical measurements.",
    icon: "science",
    to: "/proteins",
  },
  {
    title: "Mutation analysis",
    description: "Compare biological sequences and inspect sequence-level changes.",
    icon: "compare_arrows",
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
.hero {
  position: relative;
  display: flex;
  min-height: 330px;
  align-items: center;
  overflow: hidden;
  padding: 52px 56px;
  border-radius: 10px;
  background: #123d59;
  color: white;
}

.hero .biostar-eyebrow {
  color: #6ed1ca;
}

.hero .biostar-subtitle {
  color: rgb(255 255 255 / 68%);
}

.hero-copy {
  position: relative;
  z-index: 1;
  max-width: 730px;
}

.hero .biostar-title {
  color: white;
}

.hero-art {
  position: absolute;
  top: 0;
  right: 0;
  width: 38%;
  height: 100%;
  opacity: 0.55;
}

.helix {
  position: absolute;
  top: 42px;
  right: 80px;
  width: 130px;
  height: 250px;
  transform: rotate(16deg);
}

.helix::before,
.helix::after {
  position: absolute;
  top: 0;
  bottom: 0;
  width: 3px;
  content: "";
  background: #55c7c0;
  border-radius: 4px;
}

.helix::before { left: 20px; transform: rotate(8deg); }
.helix::after { right: 20px; transform: rotate(-8deg); }

.helix i {
  position: absolute;
  top: calc((var(--n) - 1) * 38px + 8px);
  left: 30px;
  width: 70px;
  height: 2px;
  background: rgb(255 255 255 / 55%);
  transform: rotate(calc((var(--n) - 4) * 7deg));
}

.status-line {
  display: flex;
  align-items: center;
  gap: 9px;
  margin: 20px 2px 54px;
  color: var(--bio-muted);
  font-size: 0.78rem;
}

.status-marker {
  width: 7px;
  height: 7px;
  border-radius: 50%;
  background: #d6a23d;
}

.status-marker.connected { background: #55a863; }
.status-marker.unavailable { background: #c65a56; }

.status-label {
  color: var(--bio-ink);
  font-weight: 750;
}

.tools-heading h2 {
  margin: 6px 0 20px;
  color: var(--bio-ink);
  font-size: 1.65rem;
  letter-spacing: -0.025em;
}

.tool-list {
  border-top: 1px solid var(--bio-line);
}

.tool-row {
  min-height: 98px;
  padding: 12px 8px;
  border-bottom: 1px solid var(--bio-line);
  border-radius: 0;
}

.tool-row:hover {
  background: rgb(23 107 135 / 4%);
}

.tool-index {
  width: 44px;
  color: #9aabb2;
  font-family: "Roboto Mono", "Courier New", monospace;
  font-size: 0.72rem;
}

.tool-title {
  color: var(--bio-ink);
  font-size: 1rem;
  font-weight: 700;
}

.tool-row :deep(.q-item__label--caption) {
  margin-top: 4px;
  color: var(--bio-muted);
}
</style>
