<template>
  <q-page class="biostar-page">
    <main class="biostar-content">
      <header>
        <div class="biostar-eyebrow">Analysis / Proteins</div>
        <h1 class="biostar-title">Protein analysis</h1>
        <p class="biostar-subtitle">
          Calculate physicochemical properties and composition from a protein sequence.
        </p>
      </header>

      <section class="biostar-panel q-mt-xl">
        <div class="biostar-panel__header">
          <div class="biostar-section-label">Input sequence</div>
          <div class="text-body2 text-grey-6 q-mt-xs">Protein sequences may contain whitespace and FASTA formatting.</div>
        </div>

        <div class="biostar-panel__body">
          <q-input
            v-model="sequence"
            type="textarea"
            outlined
            autogrow
            label="Protein sequence"
            hint="Whitespace is ignored · Maximum 10,000 characters"
            :disable="loading"
            input-class="biostar-sequence"
            @keydown.ctrl.enter.prevent="analyze"
          />

          <div class="row items-center q-mt-lg">
            <q-toggle v-model="fullAnalysis" label="Run complete analysis" color="primary" />
          </div>

          <div v-if="!fullAnalysis" class="q-mt-lg">
            <div class="biostar-section-label q-mb-sm">Requested properties</div>
            <div class="options-grid">
              <q-checkbox v-model="options.aminoacidsCount" label="Amino-acid count" />
              <q-checkbox v-model="options.isoelectricPoint" label="Isoelectric point" />
              <q-checkbox v-model="options.chargeAtPh" label="Charge at pH" />
              <q-checkbox v-model="options.aromaticity" label="Aromaticity" />
              <q-checkbox v-model="options.secondaryStructure" label="Secondary structure" />
              <q-checkbox v-model="options.molecularWeight" label="Molecular weight" />
              <q-checkbox v-model="options.hydrophobicIndex" label="Hydrophobic index" />
              <q-checkbox v-model="options.compositionRatio" label="Composition ratio" />
              <q-checkbox v-model="options.extinctionCoefficient" label="Extinction coefficient" />
            </div>

            <q-input
              v-if="options.chargeAtPh"
              v-model.number="chargePh"
              type="number"
              outlined
              dense
              label="pH"
              min="0"
              max="14"
              step="0.1"
              class="ph-input q-mt-md"
            />
          </div>

          <div class="row items-center justify-between q-mt-lg">
            <div class="text-caption text-grey-6">Ctrl + Enter to analyze</div>
            <q-btn
              unelevated
              color="primary"
              icon="science"
              label="Run analysis"
              :loading="loading"
              :disable="!sequence.trim()"
              @click="analyze"
            />
          </div>
        </div>
      </section>

      <q-banner v-if="error" rounded class="q-mt-md bg-negative text-white">
        <template #avatar><q-icon name="error_outline" /></template>
        {{ error }}
      </q-banner>

      <section v-if="result" class="biostar-panel q-mt-md">
        <div class="biostar-panel__header row items-center justify-between">
          <div>
            <div class="biostar-section-label">Analysis results</div>
            <div class="text-body2 text-grey-6 q-mt-xs">{{ result.length }} residues</div>
          </div>
          <q-btn flat round icon="content_copy" aria-label="Copy sequence" @click="copySequence">
            <q-tooltip>Copy sequence</q-tooltip>
          </q-btn>
        </div>

        <div class="biostar-panel__body">
          <div class="result-sequence biostar-sequence">{{ result.sequence }}</div>

          <div v-if="scalarResults.length > 0" class="result-grid q-mt-xl">
            <div v-for="item in scalarResults" :key="item.label" class="metric">
              <div class="metric__label">{{ item.label }}</div>
              <div class="metric__value">{{ item.value }}</div>
            </div>
          </div>

          <div v-if="result.charge_at_pH" class="result-section">
            <div class="biostar-section-label">Charge at pH {{ result.charge_at_pH.pH }}</div>
            <div class="result-value">{{ formatValue(result.charge_at_pH.charge) }}</div>
          </div>

          <div v-if="result.aminoacids_count" class="result-section">
            <div class="biostar-section-label">Amino-acid count</div>
            <div class="chip-grid">
              <q-chip v-for="(count, aminoAcid) in result.aminoacids_count" :key="aminoAcid" square>
                <span class="biostar-sequence">{{ aminoAcid }}</span>&nbsp; {{ count }}
              </q-chip>
            </div>
          </div>

          <div v-if="result.composition_ratio" class="result-section">
            <div class="biostar-section-label">Composition ratio</div>
            <div class="chip-grid">
              <q-chip v-for="(ratio, aminoAcid) in result.composition_ratio" :key="aminoAcid" square>
                <span class="biostar-sequence">{{ aminoAcid }}</span>&nbsp; {{ formatValue(ratio) }}
              </q-chip>
            </div>
          </div>

          <div v-if="result.secondary_structure_propensity" class="result-section">
            <div class="biostar-section-label">Secondary-structure propensity</div>
            <pre>{{ formatValue(result.secondary_structure_propensity) }}</pre>
          </div>

          <div v-if="result.extinction_coefficient" class="result-section">
            <div class="biostar-section-label">Extinction coefficient</div>
            <pre>{{ formatValue(result.extinction_coefficient) }}</pre>
          </div>
        </div>
      </section>
    </main>
  </q-page>
</template>

<script setup lang="ts">
import { computed, reactive, ref } from "vue";
import { Notify } from "quasar";
import { analyzeProtein } from "../services/api";
import type { ProteinAnalysisResponse, ProteinAnalysisRequest } from "../services/api";

const sequence = ref<string>("");
const loading = ref<boolean>(false);
const error = ref<string>("");
const result = ref<ProteinAnalysisResponse | null>(null);
const fullAnalysis = ref<boolean>(true);
const chargePh = ref<number>(7);

const options = reactive({
  aminoacidsCount: false,
  isoelectricPoint: false,
  chargeAtPh: false,
  aromaticity: false,
  secondaryStructure: false,
  molecularWeight: false,
  hydrophobicIndex: false,
  compositionRatio: false,
  extinctionCoefficient: false,
});

const scalarResults = computed(() => {
  if (result.value === null) return [];

  const items: Array<{ label: string; value: string }> = [];
  if (result.value.isoelectric_point !== undefined) items.push({ label: "Isoelectric point", value: formatValue(result.value.isoelectric_point) });
  if (result.value.aromaticity !== undefined) items.push({ label: "Aromaticity", value: formatValue(result.value.aromaticity) });
  if (result.value.molecular_weight !== undefined) items.push({ label: "Molecular weight", value: formatValue(result.value.molecular_weight) });
  if (result.value.hydrophobic_index !== undefined) items.push({ label: "Hydrophobic index", value: formatValue(result.value.hydrophobic_index) });
  return items;
});

const buildRequest = (): ProteinAnalysisRequest => ({
  sequence: sequence.value,
  get_full_test_results: fullAnalysis.value,
  get_aminoacids_count: options.aminoacidsCount,
  get_isoelectric_point: options.isoelectricPoint,
  get_charge_at_pH: !fullAnalysis.value && options.chargeAtPh ? chargePh.value : null,
  get_aromaticity: options.aromaticity,
  get_secondary_structure_propensity: options.secondaryStructure,
  get_molecular_weight: options.molecularWeight,
  get_hydrophobic_index: options.hydrophobicIndex,
  get_composition_ratio: options.compositionRatio,
  get_extinction_coefficient: options.extinctionCoefficient,
});

const analyze = async (): Promise<void> => {
  if (!sequence.value.trim()) return;

  loading.value = true;
  error.value = "";
  result.value = null;

  try {
    const response = await analyzeProtein(buildRequest());

    if (response.data === null) {
      error.value = response.message?.message ?? "Protein analysis returned no data.";
      return;
    }

    result.value = response.data;
  } catch (requestError: unknown) {
    error.value = requestError instanceof Error ? requestError.message : "Protein analysis failed.";
  } finally {
    loading.value = false;
  }
};

const copySequence = async (): Promise<void> => {
  if (!result.value) return;
  try {
    await navigator.clipboard.writeText(result.value.sequence);
    Notify.create({ type: "positive", message: "Protein sequence copied." });
  } catch {
    Notify.create({ type: "negative", message: "Unable to copy the sequence." });
  }
};

const formatValue = (value: unknown): string => {
  if (typeof value === "number") return Number.isInteger(value) ? String(value) : value.toFixed(4);
  if (typeof value === "string") return value;
  return JSON.stringify(value) ?? String(value);
};
</script>

<style scoped>
.options-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 8px 16px;
}

.ph-input {
  max-width: 180px;
}

.result-sequence {
  padding: 16px;
  overflow-x: auto;
  border: 1px solid var(--biostar-border);
  border-radius: 4px;
  background: var(--biostar-bg);
  line-height: 1.7;
  word-break: break-word;
}

.result-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 12px;
}

.metric {
  padding: 18px;
  border: 1px solid var(--biostar-border);
  border-radius: 4px;
}

.metric__label {
  color: var(--biostar-muted);
  font-size: 0.75rem;
  font-weight: 700;
  letter-spacing: 0.04em;
  text-transform: uppercase;
}

.metric__value {
  margin-top: 6px;
  color: var(--biostar-text);
  font-size: 1.25rem;
  font-weight: 700;
}

.result-section {
  margin-top: 28px;
  padding-top: 24px;
  border-top: 1px solid var(--biostar-border);
}

.result-value {
  margin-top: 8px;
  font-size: 1.1rem;
}

.chip-grid {
  display: flex;
  flex-wrap: wrap;
  gap: 4px;
  margin-top: 10px;
}

pre {
  margin: 10px 0 0;
  padding: 14px;
  overflow-x: auto;
  border: 1px solid var(--biostar-border);
  background: var(--biostar-bg);
  white-space: pre-wrap;
}
</style>
