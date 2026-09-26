<template>
  <q-layout view="hHh Lpr lFf">
    <q-header class="site-header">
      <q-toolbar class="site-toolbar">
        <q-btn flat round dense icon="menu" class="menu-button" @click="drawerOpen = !drawerOpen" />
        <div class="brand">
          <div class="brand-mark">
            <span></span><span></span><span></span>
          </div>
          <div>
            <div class="brand-name">BioSTAR</div>
            <div class="brand-caption">BIOINFORMATICS RESOURCE</div>
          </div>
        </div>
        <q-space />
        <q-btn flat round icon="dark_mode" class="theme-button" @click="toggleDarkMode" />
      </q-toolbar>
    </q-header>

    <q-drawer v-model="drawerOpen" show-if-above :width="232" class="site-drawer">
      <div class="drawer-inner">
        <div class="drawer-title">TOOLS</div>

        <q-item clickable v-ripple to="/" exact class="nav-link">
          <q-item-section avatar><q-icon name="home" /></q-item-section>
          <q-item-section>Overview</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/sequences" class="nav-link">
          <q-item-section avatar><q-icon name="biotech" /></q-item-section>
          <q-item-section>Sequences</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/proteins" class="nav-link">
          <q-item-section avatar><q-icon name="science" /></q-item-section>
          <q-item-section>Proteins</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/mutations" class="nav-link">
          <q-item-section avatar><q-icon name="compare_arrows" /></q-item-section>
          <q-item-section>Mutations</q-item-section>
        </q-item>

        <div class="drawer-divider"></div>
        <div class="drawer-title">PLATFORM</div>

        <q-item clickable v-ripple to="/settings" class="nav-link">
          <q-item-section avatar><q-icon name="settings" /></q-item-section>
          <q-item-section>Settings</q-item-section>
        </q-item>

        <q-space />
        <div class="drawer-footer">
          <div>BioSTAR API</div>
          <span>Scientific computing tools</span>
        </div>
      </div>
    </q-drawer>

    <q-page-container>
      <router-view />
    </q-page-container>
  </q-layout>
</template>

<script setup lang="ts">
import { ref } from "vue";
import { Dark } from "quasar";

const drawerOpen = ref<boolean>(true);

const toggleDarkMode = (): void => {
  Dark.toggle();
};
</script>

<style scoped>
.site-header {
  background: linear-gradient(115deg, #123d59 0%, #176b87 58%, #238b8f 100%);
  color: white;
}

.site-toolbar {
  min-height: 4.25rem;
  padding: 0 2%;
}

.menu-button,
.theme-button {
  color: rgb(255 255 255 / 78%);
}

.brand {
  display: flex;
  align-items: center;
  gap: 0.7rem;
  margin-left: 0.75rem;
}

.brand-mark {
  display: flex;
  width: 1.875rem;
  height: 1.875rem;
  align-items: flex-end;
  gap: 0.125rem;
}

.brand-mark span {
  width: 0.4375rem;
  border-radius: 0.3rem 0.3rem 0.125rem 0.125rem;
  background: #5bc6c0;
  transform: skewY(-22deg);
}

.brand-mark span:nth-child(1) { height: 0.95rem; }
.brand-mark span:nth-child(2) { height: 1.45rem; }
.brand-mark span:nth-child(3) { height: 1.875rem; }

.brand-name {
  font-size: 1.15rem;
  font-weight: 800;
  letter-spacing: 0.06em;
}

.brand-caption {
  margin-top: 0.0625rem;
  color: rgb(255 255 255 / 58%);
  font-size: 0.56rem;
  letter-spacing: 0.13em;
}

.site-drawer {
  background: #fff;
  border-right: 1px solid #e4eaee;
}

.drawer-inner {
  display: flex;
  min-height: 100%;
  flex-direction: column;
  padding: 1.75rem 0.75rem 1.125rem;
}

.drawer-title {
  padding: 0 0.8rem 0.5rem;
  color: #82929c;
  font-size: 0.65rem;
  font-weight: 800;
  letter-spacing: 0.14em;
}

.nav-link {
  min-height: 2.75rem;
  margin: 0.125rem 0;
  border-radius: 0.35rem;
  color: #435963;
  font-size: 0.9rem;
}

.nav-link .q-icon {
  color: #82949d;
}

.nav-link.q-router-link--active {
  background: linear-gradient(100deg, #e3f3f4 0%, #edf8f8 100%);
  color: #126d79;
  font-weight: 700;
}

.nav-link.q-router-link--active .q-icon {
  color: #168391;
}

.drawer-divider {
  height: 1px;
  margin: 1.125rem 0.75rem;
  background: #e4eaee;
}

.drawer-footer {
  margin: 1.125rem 0.75rem 0;
  padding-top: 0.875rem;
  border-top: 1px solid #e4eaee;
  color: #536872;
  font-size: 0.75rem;
  font-weight: 700;
}

.drawer-footer span {
  display: block;
  margin-top: 0.1875rem;
  color: #91a0a8;
  font-size: 0.68rem;
  font-weight: 400;
}

@media (max-width: 37.49rem) {
  .site-toolbar {
    min-height: 3.75rem;
  }

  .brand-caption {
    display: none;
  }

  .brand {
    gap: 0.55rem;
    margin-left: 0.35rem;
  }

  .brand-name {
    font-size: 1rem;
  }
}
