<template>
  <q-page class="biostar-page">
    <main class="biostar-content">
      <section class="welcome">
        <div class="welcome-copy">
          <div class="biostar-eyebrow">BioSTAR workspace</div>
          <h1 class="biostar-title">What are you<br class="gt-xs" ></br> working on?</h1>
          <p class="biostar-subtitle">
            Run focused sequence and protein analyses through one scientific workspace.
          </p>
        </div>

        <div class="welcome-orbit" aria-hidden="true">
          <div class="orbit orbit--outer"></div>
          <div class="orbit orbit--inner"></div>
          <span class="orbit-dot orbit-dot--one"></span>
          <span class="orbit-dot orbit-dot--two"></span>
          <span class="orbit-dot orbit-dot--three"></span>
        </div>
      </section>

      <section class="workspace-meta">
        <div class="api-state">
          <span class="state-dot" :class="apiStatus"></span>
          <span class="api-name">API</span>
          <span>{{ apiStatusLabel }}</span>
        </div>
        <span class="meta-separator">·</span>
        <span>Scientific computing engine</span>
      </section>

      <section class="tools-section">
        <header class="section-header">
          <div>
            <div class="biostar-eyebrow">Tools</div>
            <h2>Start an analysis</h2>
          </div>
          <span class="tool-count">{{ cards.length }} workflows</span>
        </header>

        <div class="tool-grid">
          <q-card
            v-for="card in cards"
            :key="card.title"
            flat
            bordered
            class="tool-card"
            clickable
            @click="$router.push(card.to)"
          >
            <q-card-section class="tool-card__top">
              <div class="tool-icon">
                <q-icon :name="card.icon" ></q-icon>
              </div>
              <q-icon name="arrow_outward" class="tool-arrow" ></q-icon>
            </q-card-section>
            <q-card-section class="tool-card__body">
              <div class="tool-number">{{ card.number }}</div>
              <h3>{{ card.title }}</h3>
              <p>{{ card.description }}</p>
            </q-card-section>
            <q-card-section class="tool-card__footer">
              <span>Open workspace</span>
              <q-icon name="arrow_forward" ></q-icon>
            </q-card-section>
          </q-card>
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
  number: string;
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
    number: "01",
    title: "Sequence conversion",
    description: "Convert DNA and RNA sequences, or translate nucleotides into proteins.",
    icon: "biotech",
    to: "/sequences",
  },
  {
    number: "02",
    title: "Protein analysis",
    description: "Explore physicochemical properties, composition and biochemical measurements.",
    icon: "science",
    to: "/proteins",
  },
  {
    number: "03",
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
.welcome {
  position: relative;
  display: flex;
  min-height: 21rem;
  align-items: center;
  overflow: hidden;
  padding: clamp(2rem, 5vw, 4rem);
  border: 1px solid rgb(20 125 131 / 18%);
  border-radius: 1.25rem;
  background:
    radial-gradient(circle at 88% 40%, rgb(83 184 174 / 22%), transparent 18rem),
    linear-gradient(135deg, #103b4b 0%, #126b73 58%, #167f7e 100%);
  box-shadow: var(--bio-shadow);
  color: white;
}

.title-break {
  display: block;
}

.welcome-copy {
  position: relative;
  z-index: 2;
  max-width: 42rem;
}

.welcome .biostar-eyebrow {
  color: #7bd2ca;
}

.welcome .biostar-title {
  color: white;
}

.welcome .biostar-subtitle {
  max-width: 38rem;
  color: rgb(255 255 255 / 72%);
}

.welcome-orbit {
  position: absolute;
  top: 50%;
  right: 8%;
  width: clamp(10rem, 25vw, 18rem);
  aspect-ratio: 1;
  transform: translateY(-50%);
}

.orbit {
  position: absolute;
  inset: 0;
  border: 1px solid rgb(255 255 255 / 20%);
  border-radius: 50%;
}

.orbit--inner {
  inset: 16%;
  border-color: rgb(123 210 202 / 34%);
  transform: rotate(28deg) scaleX(0.48);
}

.orbit--outer {
  transform: rotate(-28deg) scaleX(0.55);
}

.orbit-dot {
  position: absolute;
  width: 0.65rem;
  height: 0.65rem;
  border-radius: 50%;
  background: #78d5cc;
  box-shadow: 0 0 0 0.4rem rgb(120 213 204 / 10%);
}

.orbit-dot--one { top: 15%; left: 28%; }
.orbit-dot--two { right: 5%; bottom: 24%; }
.orbit-dot--three { bottom: 8%; left: 44%; }

.workspace-meta {
  display: flex;
  align-items: center;
  gap: 0.55rem;
  margin: 1rem 0 3.25rem;
  color: var(--bio-muted);
  font-size: 0.75rem;
}

.api-state {
  display: flex;
  align-items: center;
  gap: 0.45rem;
}

.state-dot {
  width: 0.45rem;
  height: 0.45rem;
  border-radius: 50%;
  background: #d4a33e;
}

.state-dot.connected { background: #49a36f; }
.state-dot.unavailable { background: #c75c59; }

.api-name {
  color: var(--bio-ink);
  font-weight: 750;
}

.meta-separator {
  color: var(--bio-line);
}

.section-header {
  display: flex;
  align-items: flex-end;
  justify-content: space-between;
  gap: 1rem;
  margin-bottom: 1.25rem;
}

.section-header h2 {
  margin: 0.4rem 0 0;
  color: var(--bio-ink);
  font-size: 1.55rem;
  letter-spacing: -0.025em;
}

.tool-count {
  color: var(--bio-muted);
  font-size: 0.72rem;
}

.tool-grid {
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: 1rem;
}

.tool-card {
  min-height: 18rem;
  overflow: hidden;
  border-color: var(--bio-line);
  border-radius: 1rem;
  background: var(--bio-paper);
  box-shadow: none;
  transition: transform 160ms ease, border-color 160ms ease, box-shadow 160ms ease;
}

.tool-card:hover {
  transform: translateY(-0.2rem);
  border-color: rgb(20 125 131 / 35%);
  box-shadow: var(--bio-shadow-small);
}

.tool-card__top {
  display: flex;
  align-items: flex-start;
  justify-content: space-between;
  padding: 1.25rem 1.25rem 0;
}

.tool-icon {
  display: grid;
  width: 2.8rem;
  height: 2.8rem;
  place-items: center;
  border-radius: 0.8rem;
  background: var(--bio-soft);
  color: var(--bio-primary);
  font-size: 1.35rem;
}

.tool-arrow {
  color: var(--bio-muted);
}

.tool-card__body {
  padding: 1rem 1.25rem 1.25rem;
}

.tool-number {
  color: var(--bio-muted);
  font-family: "Roboto Mono", "Courier New", monospace;
  font-size: 0.68rem;
}

.tool-card h3 {
  margin: 0.5rem 0 0;
  color: var(--bio-ink);
  font-size: 1.08rem;
  font-weight: 750;
}

.tool-card p {
  margin: 0.65rem 0 0;
  color: var(--bio-muted);
  font-size: 0.86rem;
  line-height: 1.6;
}

.tool-card__footer {
  display: flex;
  align-items: center;
  justify-content: space-between;
  margin-top: auto;
  padding: 0.9rem 1.25rem;
  border-top: 1px solid var(--bio-line);
  color: var(--bio-primary);
  font-size: 0.72rem;
  font-weight: 700;
}

@media (max-width: 59.99rem) {
  .welcome-orbit {
    right: -2%;
  }

  .tool-grid {
    grid-template-columns: 1fr;
  }

  .tool-card {
    min-height: 14rem;
  }
}

@media (max-width: 37.49rem) {
  .title-break {
    display: inline;
  }

  .welcome {
    min-height: 23rem;
    align-items: flex-start;
    padding: 2rem 1.25rem;
  }

  .welcome-orbit {
    top: auto;
    right: -8%;
    bottom: -18%;
    width: 12rem;
    transform: none;
    opacity: 0.8;
  }

  .workspace-meta {
    margin-bottom: 2.5rem;
  }

  .section-header {
    align-items: flex-start;
  }
}
