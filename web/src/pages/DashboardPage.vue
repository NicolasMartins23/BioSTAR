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
  min-height: 20rem;
  align-items: center;
  overflow: hidden;
  padding: 3.25rem 5%;
  border-radius: 0.9rem;
  background: linear-gradient(115deg, #123d59 0%, #176b87 62%, #238b8f 100%);
  box-shadow: 0 1.5rem 3rem rgb(18 61 89 / 14%);
  color: white;
}

.hero .biostar-eyebrow {
  color: #6ed1ca;
}

.hero .biostar-subtitle {
  color: rgb(255 255 255 / 70%);
}

.hero-copy {
  position: relative;
  z-index: 1;
  max-width: 45rem;
}

.hero .biostar-title {
  color: white;
}

.hero-art {
  position: absolute;
  inset: 0 0 0 auto;
  width: 38%;
  opacity: 0.5;
}

.helix {
  position: absolute;
  top: 12%;
  right: 16%;
  width: 8rem;
  height: 15rem;
  transform: rotate(16deg);
}

.helix::before,
.helix::after {
  position: absolute;
  top: 0;
  bottom: 0;
  width: 0.2rem;
  content: "";
  background: #55c7c0;
  border-radius: 0.25rem;
}

.helix::before { left: 12%; transform: rotate(8deg); }
.helix::after { right: 12%; transform: rotate(-8deg); }

.helix i {
  position: absolute;
  top: calc((var(--n) - 1) * 2.375rem + 0.5rem);
  left: 23%;
  width: 54%;
  height: 0.125rem;
  background: rgb(255 255 255 / 55%);
  transform: rotate(calc((var(--n) - 4) * 7deg));
}

.status-line {
  display: flex;
  align-items: center;
  gap: 0.55rem;
  margin: 1.25rem 0.125rem 3.5rem;
  color: var(--bio-muted);
  font-size: 0.78rem;
}

.status-marker {
  width: 0.45rem;
  height: 0.45rem;
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
  margin: 0.375rem 0 1.25rem;
  color: var(--bio-ink);
  font-size: 1.65rem;
  letter-spacing: -0.025em;
}

.tool-list {
  border-top: 1px solid var(--bio-line);
}

.tool-row {
  min-height: 6rem;
  padding: 0.75rem 0.5rem;
  border-bottom: 1px solid var(--bio-line);
  border-radius: 0;
  transition: background 160ms ease, padding 160ms ease;
}

.tool-row:hover {
  background: linear-gradient(90deg, rgb(23 107 135 / 6%), rgb(40 165 160 / 2%));
  padding-inline: 0.75rem;
}

.tool-index {
  width: 2.75rem;
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
  margin-top: 0.25rem;
  color: var(--bio-muted);
}

@media (max-width: 59.99rem) {
  .hero {
    min-height: 18rem;
    padding: 2.75rem 5%;
  }

  .hero-art {
    width: 30%;
  }

  .status-line {
    margin-bottom: 2.75rem;
  }
}

@media (max-width: 37.49rem) {
  .hero {
    min-height: 25rem;
    align-items: flex-start;
    padding: 2rem 1.25rem;
  }

  .hero .biostar-title {
    font-size: clamp(2rem, 10vw, 2.8rem);
  }

  .hero .biostar-subtitle {
    max-width: 100%;
    font-size: 0.92rem;
  }

  .hero-art {
    top: auto;
    right: -10%;
    bottom: -20%;
    width: 75%;
    height: 60%;
    opacity: 0.28;
  }

  .helix {
    top: 0;
    right: 12%;
  }

  .status-line {
    margin: 1rem 0.125rem 2.5rem;
  }

  .tool-row {
    min-height: 5.5rem;
    padding: 0.75rem 0;
  }

  .tool-row:hover {
    padding-inline: 0.25rem;
  }

  .tool-index {
    display: none;
  }

  .tool-row :deep(.q-item__section--avatar) {
    min-width: 2.75rem;
  }

  .tool-row :deep(.q-item__section--side:last-child) {
    padding-left: 0.5rem;
  }
}

</style>
