<template>
  <q-page class="biostar-page">
    <main class="biostar-content">
      <header>
        <div class="biostar-eyebrow">Analysis / Nucleic acids</div>
        <h1 class="biostar-title">Sequence conversion</h1>
        <p class="biostar-subtitle">
          Convert DNA and RNA sequences using the BioSTAR biochemical engine.
          Single-sequence FASTA input is supported.
        </p>
      </header>

      <section class="biostar-panel q-mt-xl">
        <div class="biostar-panel__header">
          <div class="biostar-section-label">Input</div>
          <div class="text-body2 text-grey-6 q-mt-xs">Define the conversion and provide a sequence.</div>
        </div>
        <div class="biostar-panel__body">
          <div class="row q-col-gutter-md">
            <div class="col-12 col-md-4">
              <q-select
                v-model="conversion"
                :options="conversionOptions"
                emit-value
                map-options
                outlined
                label="Conversion"
                :disable="isLoading"
              />
            </div>
            <div class="col-12 col-md-8">
              <q-input
                v-model="sequence"
                outlined
                type="textarea"
                autogrow
                label="Sequence"
                hint="Whitespace is ignored · Maximum 1,000 characters"
                :disable="isLoading"
                input-class="biostar-sequence"
                @keydown.ctrl.enter="convert"
              />
            </div>
          </div>

          <div class="row items-center justify-between q-mt-lg">
            <div class="text-caption text-grey-6">Ctrl + Enter to convert</div>
            <q-btn
              unelevated
              color="primary"
              icon="play_arrow"
              label="Run conversion"
              :loading="isLoading"
              :disable="!sequence.trim()"
              @click="convert"
            />
          </div>
        </div>
      </section>

      <q-banner v-if="errorMessage" rounded class="q-mt-md bg-negative text-white">
        <template #avatar><q-icon name="error_outline" /></template>
        {{ errorMessage }}
      </q-banner>

      <section v-if="result" class="biostar-panel q-mt-md">
        <div class="biostar-panel__header row items-center justify-between">
          <div>
            <div class="biostar-section-label">Result</div>
            <div class="text-body2 text-grey-6 q-mt-xs">{{ conversionLabel }}</div>
          </div>
          <q-btn flat round icon="content_copy" aria-label="Copy result" @click="copyResult">
            <q-tooltip>Copy result</q-tooltip>
          </q-btn>
        </div>
        <div class="biostar-panel__body">
          <div class="result-box biostar-sequence">{{ result }}</div>
        </div>
      </section>
    </main>
  </q-page>
</template>

<script setup lang="ts">
import { computed, ref } from "vue";
import { Notify } from "quasar";
import type { APIResponse, SequenceResponse } from "../services/api";
import {
  convertDnaToProtein,
  convertDnaToRna,
  convertRnaToDna,
  convertRnaToProtein,
} from "../services/api";

type Conversion = "dna-rna" | "dna-protein" | "rna-protein" | "rna-dna";

const conversion = ref<Conversion>("dna-rna");
const sequence = ref<string>("");
const result = ref<string>("");
const errorMessage = ref<string>("");
const isLoading = ref<boolean>(false);

const conversionOptions: Array<{ label: string; value: Conversion }> = [
  { label: "DNA → RNA", value: "dna-rna" },
  { label: "DNA → Protein", value: "dna-protein" },
  { label: "RNA → Protein", value: "rna-protein" },
  { label: "RNA → DNA", value: "rna-dna" },
];

const conversionLabel = computed<string>(() => {
  const option = conversionOptions.find((item) => item.value === conversion.value);
  return option?.label ?? conversion.value;
});

const convert = async (): Promise<void> => {
  const input: string = sequence.value.trim();
  if (!input) return;

  isLoading.value = true;
  errorMessage.value = "";
  result.value = "";

  try {
    let response: APIResponse<SequenceResponse>;

    if (conversion.value === "dna-rna") response = await convertDnaToRna(input);
    else if (conversion.value === "dna-protein") response = await convertDnaToProtein(input);
    else if (conversion.value === "rna-protein") response = await convertRnaToProtein(input);
    else response = await convertRnaToDna(input);

    if (response.data === null) {
      errorMessage.value = response.message?.message ?? "Sequence conversion returned no data.";
      return;
    }

    result.value = response.data.sequence;
  } catch (error: unknown) {
    errorMessage.value = error instanceof Error ? error.message : "Sequence conversion failed.";
  } finally {
    isLoading.value = false;
  }
};

const copyResult = async (): Promise<void> => {
  if (!result.value) return;

  try {
    await navigator.clipboard.writeText(result.value);
    Notify.create({ type: "positive", message: "Result copied to clipboard." });
  } catch {
    Notify.create({ type: "negative", message: "Unable to copy the result." });
  }
};
</script>

<style scoped>
.result-box {
  min-height: 110px;
  padding: 18px;
  overflow-x: auto;
  border: 1px solid var(--biostar-border);
  border-radius: 4px;
  background: var(--biostar-bg);
  white-space: pre-wrap;
  word-break: break-word;
  line-height: 1.7;
}
</style>
