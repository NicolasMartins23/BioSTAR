<template>
  <q-page padding>
    <div class="page-shell">
      <div class="text-overline text-primary">Analysis</div>
      <div class="text-h4 text-weight-bold">Proteins</div>
      <div class="text-body1 text-grey-6 q-mt-sm">
        Analyze a protein sequence with the BioSTAR biochemical engine.
      </div>

      <q-card flat bordered class="q-mt-xl">
        <q-card-section>
          <q-input v-model="sequence" type="textarea" outlined autogrow
            label="Protein sequence"
            hint="Whitespace is ignored. Maximum 10,000 characters."
            :disable="loading"
            @keydown.ctrl.enter.prevent="analyze" />

          <q-toggle v-model="fullAnalysis" label="Full analysis" color="primary" class="q-mt-md" />

          <div v-if="!fullAnalysis" class="options-grid q-mt-md">
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

          <q-input v-if="!fullAnalysis && options.chargeAtPh"
            v-model.number="chargePh" type="number" outlined dense
            label="pH" min="0" max="14" step="0.1" class="ph-input q-mt-md" />

          <q-btn color="primary" icon="science" label="Analyze protein"
            :loading="loading" :disable="!sequence.trim()"
            class="q-mt-lg" @click="analyze" />
        </q-card-section>
      </q-card>

      <q-banner v-if="error" rounded class="bg-negative text-white q-mt-lg">
        {{ error }}
      </q-banner>

      <q-card v-if="result" flat bordered class="q-mt-lg">
        <q-card-section>
          <div class="row items-center justify-between">
            <div>
              <div class="text-h6">Results</div>
              <div class="text-caption text-grey-6">Length: {{ result.length }}</div>
            </div>
            <q-btn flat round icon="content_copy" @click="copySequence" />
          </div>

          <q-input :model-value="result.sequence" readonly outlined type="textarea"
            autogrow class="q-mt-md" />
        </q-card-section>

        <q-separator />

        <q-card-section v-if="scalarResults.length > 0">
          <div class="result-grid">
            <q-card v-for="item in scalarResults" :key="item.label" flat bordered>
              <q-card-section>
                <div class="text-caption text-grey-6">{{ item.label }}</div>
                <div class="text-h6 q-mt-xs">{{ item.value }}</div>
              </q-card-section>
            </q-card>
          </div>
        </q-card-section>

        <q-card-section v-if="result.charge_at_pH">
          <div class="text-subtitle1">Charge at pH {{ result.charge_at_pH.pH }}</div>
          <div class="text-body1">{{ formatValue(result.charge_at_pH.charge) }}</div>
        </q-card-section>

        <q-card-section v-if="result.aminoacids_count">
          <div class="text-subtitle1 q-mb-sm">Amino-acid count</div>
          <q-chip v-for="(count, aminoAcid) in result.aminoacids_count"
            :key="aminoAcid" square>
            {{ aminoAcid }}: {{ count }}
          </q-chip>
        </q-card-section>

        <q-card-section v-if="result.composition_ratio">
          <div class="text-subtitle1 q-mb-sm">Composition ratio</div>
          <q-chip v-for="(ratio, aminoAcid) in result.composition_ratio"
            :key="aminoAcid" square>
            {{ aminoAcid }}: {{ formatValue(ratio) }}
          </q-chip>
        </q-card-section>

        <q-card-section v-if="result.secondary_structure_propensity">
          <div class="text-subtitle1">Secondary-structure propensity</div>
          <pre>{{ formatValue(result.secondary_structure_propensity) }}</pre>
        </q-card-section>

        <q-card-section v-if="result.extinction_coefficient">
          <div class="text-subtitle1">Extinction coefficient</div>
          <pre>{{ formatValue(result.extinction_coefficient) }}</pre>
        </q-card-section>
      </q-card>
    </div>
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
  await navigator.clipboard.writeText(result.value.sequence);
  Notify.create({ type: "positive", message: "Protein sequence copied." });
};

const formatValue = (value: unknown): string => {
  if (typeof value === "number") return Number.isInteger(value) ? String(value) : value.toFixed(4);
  if (typeof value === "string") return value;
  return JSON.stringify(value);
};
</script>

<style scoped>
.page-shell { max-width: 1200px; margin: 0 auto; }
.options-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 8px 16px;
}
.ph-input { max-width: 180px; }
.result-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 12px;
}
pre { white-space: pre-wrap; overflow-x: auto; }
</style>
