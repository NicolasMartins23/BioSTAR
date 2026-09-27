<template>
  <q-layout view="hHh Lpr lFf">
    <q-header class="site-header">
      <q-toolbar class="site-toolbar">
        <q-btn flat round dense icon="menu" class="menu-button" @click="drawerOpen = !drawerOpen" ></q-btn>
        <router-link to="/" class="brand">
          <span class="brand-mark" aria-hidden="true">
            <i></i><i></i><i></i>
          </span>
          <span class="brand-copy">
            <strong>BioSTAR</strong>
            <small>Biological analysis platform</small>
          </span>
        </router-link>
        <q-space></q-space>
        <div class="header-context gt-xs">Scientific computing workspace</div>
        <q-btn flat round icon="dark_mode" class="theme-button" @click="toggleDarkMode" ></q-btn>
      </q-toolbar>
    </q-header>

    <q-drawer v-model="drawerOpen" show-if-above :width="248" class="site-drawer">
      <div class="drawer-inner">
        <div class="drawer-brand">
          <div class="drawer-brand__mark" aria-hidden="true">
            <i></i><i></i><i></i>
          </div>
          <div>
            <strong>BioSTAR</strong>
            <span>Scientific workspace</span>
          </div>
        </div>

        <div class="drawer-intro">
          <span>WORKSPACE</span>
          <strong>Analysis tools</strong>
        </div>

        <nav class="nav-group" aria-label="Analysis tools">
          <q-item clickable v-ripple to="/" exact class="nav-link">
            <q-item-section avatar><span class="nav-icon"><q-icon name="space_dashboard" ></q-icon></span></q-item-section>
            <q-item-section>Overview</q-item-section>
            <q-item-section side class="nav-active-mark"><span></span></q-item-section>
          </q-item>
          <q-item clickable v-ripple to="/sequences" class="nav-link">
            <q-item-section avatar><span class="nav-icon"><q-icon name="biotech" ></q-icon></span></q-item-section>
            <q-item-section>Sequences</q-item-section>
            <q-item-section side class="nav-active-mark"><span></span></q-item-section>
          </q-item>
          <q-item clickable v-ripple to="/proteins" class="nav-link">
            <q-item-section avatar><span class="nav-icon"><q-icon name="science" ></q-icon></span></q-item-section>
            <q-item-section>Proteins</q-item-section>
            <q-item-section side class="nav-active-mark"><span></span></q-item-section>
          </q-item>
          <q-item clickable v-ripple to="/mutations" class="nav-link">
            <q-item-section avatar><span class="nav-icon"><q-icon name="compare_arrows" ></q-icon></span></q-item-section>
            <q-item-section>Mutations</q-item-section>
            <q-item-section side class="nav-active-mark"><span></span></q-item-section>
          </q-item>
        </nav>

        <div class="drawer-divider"></div>

        <div class="drawer-intro drawer-intro--small">
          <span>PLATFORM</span>
        </div>

        <q-item clickable v-ripple to="/settings" class="nav-link">
          <q-item-section avatar><span class="nav-icon"><q-icon name="tune" ></q-icon></span></q-item-section>
          <q-item-section>Settings</q-item-section>
          <q-item-section side class="nav-active-mark"><span></span></q-item-section>
        </q-item>

        <q-space></q-space>

        <div class="drawer-status">
          <span class="status-dot"></span>
          <div>
            <strong>BioSTAR API</strong>
            <small>Scientific engine</small>
          </div>
        </div>
      </div>
    </q-drawer>

    <q-page-container>
      <router-view></router-view>
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
  background: rgb(255 255 255 / 92%);
  border-bottom: 1px solid var(--bio-line);
  color: var(--bio-ink);
  backdrop-filter: blur(1rem);
}

.body--dark .site-header {
  background: rgb(13 21 24 / 92%);
}

.site-toolbar {
  min-height: 4rem;
  padding-inline: 1.25rem;
}

.menu-button,
.theme-button {
  color: var(--bio-muted);
}

.brand {
  display: flex;
  align-items: center;
  gap: 0.7rem;
  margin-left: 0.75rem;
  color: inherit;
  text-decoration: none;
}

.brand-mark {
  display: flex;
  width: 1.8rem;
  height: 1.8rem;
  align-items: flex-end;
  gap: 0.18rem;
}

.brand-mark i {
  display: block;
  width: 0.42rem;
  border-radius: 0.4rem;
  background: var(--bio-primary);
  transform: skewY(-18deg);
}

.brand-mark i:nth-child(1) { height: 0.85rem; }
.brand-mark i:nth-child(2) { height: 1.3rem; }
.brand-mark i:nth-child(3) { height: 1.7rem; }

.brand-copy {
  display: flex;
  flex-direction: column;
  line-height: 1;
}

.brand-copy strong {
  font-size: 1.05rem;
  letter-spacing: 0.04em;
}

.brand-copy small {
  margin-top: 0.3rem;
  color: var(--bio-muted);
  font-size: 0.62rem;
  letter-spacing: 0.02em;
}

.header-context {
  margin-right: 0.75rem;
  color: var(--bio-muted);
  font-size: 0.75rem;
}

.site-drawer {
  overflow: hidden;
  background:
    radial-gradient(circle at 100% 8%, rgb(83 184 174 / 8%), transparent 13rem),
    var(--bio-paper);
  border-right: 1px solid var(--bio-line);
}

.drawer-inner {
  display: flex;
  min-height: 100%;
  flex-direction: column;
  padding: 1.25rem 0.8rem 1rem;
}

.drawer-brand {
  display: flex;
  align-items: center;
  gap: 0.7rem;
  margin: 0.25rem 0.45rem 1.75rem;
  padding-bottom: 1.25rem;
  border-bottom: 1px solid var(--bio-line);
}

.drawer-brand__mark {
  display: flex;
  width: 1.65rem;
  height: 1.65rem;
  align-items: flex-end;
  gap: 0.15rem;
}

.drawer-brand__mark i {
  display: block;
  width: 0.36rem;
  border-radius: 0.4rem;
  background: var(--bio-primary);
  transform: skewY(-18deg);
}

.drawer-brand__mark i:nth-child(1) { height: 0.72rem; }
.drawer-brand__mark i:nth-child(2) { height: 1.1rem; }
.drawer-brand__mark i:nth-child(3) { height: 1.55rem; }

.drawer-brand div:last-child {
  display: flex;
  flex-direction: column;
}

.drawer-brand strong {
  color: var(--bio-ink);
  font-size: 0.92rem;
  letter-spacing: 0.04em;
}

.drawer-brand span {
  margin-top: 0.2rem;
  color: var(--bio-muted);
  font-size: 0.62rem;
}

.drawer-intro {
  display: flex;
  flex-direction: column;
  gap: 0.35rem;
  padding: 0.35rem 0.75rem 0.7rem;
}

.drawer-intro span {
  color: var(--bio-primary);
  font-size: 0.6rem;
  font-weight: 800;
  letter-spacing: 0.14em;
}

.drawer-intro strong {
  color: var(--bio-ink);
  font-size: 0.9rem;
  font-weight: 700;
}

.drawer-intro--small {
  padding-bottom: 0.45rem;
}

.nav-group {
  position: relative;
}

.nav-link {
  position: relative;
  min-height: 3rem;
  margin: 0.18rem 0;
  border-radius: 0.7rem;
  color: var(--bio-muted);
  font-size: 0.86rem;
  transition: background 140ms ease, color 140ms ease, transform 140ms ease;
}

.nav-link .q-icon {
  color: inherit;
}

.nav-icon {
  display: grid;
  width: 2rem;
  height: 2rem;
  place-items: center;
  border-radius: 0.55rem;
  background: transparent;
  transition: background 140ms ease, color 140ms ease;
}

.nav-link:hover {
  background: var(--bio-soft);
  color: var(--bio-ink);
  transform: translateX(0.12rem);
}

.nav-link:hover .nav-icon {
  background: rgb(20 125 131 / 8%);
  color: var(--bio-primary);
}

.nav-link.q-router-link--active {
  background: linear-gradient(100deg, var(--bio-soft-strong), rgb(83 184 174 / 5%));
  color: var(--bio-primary-strong);
  font-weight: 700;
}

.nav-link.q-router-link--active .nav-icon {
  background: rgb(20 125 131 / 12%);
  color: var(--bio-primary);
}

.nav-active-mark {
  width: 0.3rem;
  min-width: 0.3rem;
  padding: 0;
}

.nav-active-mark span {
  display: block;
  width: 0.18rem;
  height: 1.35rem;
  border-radius: 0.2rem;
  background: var(--bio-primary);
  opacity: 0;
  transform: scaleY(0.5);
  transition: opacity 140ms ease, transform 140ms ease;
}

.nav-link.q-router-link--active .nav-active-mark span {
  opacity: 1;
  transform: scaleY(1);
}

.drawer-divider {
  height: 1px;
  margin: 1.15rem 0.75rem;
  background: var(--bio-line);
}

.drawer-status {
  display: flex;
  align-items: center;
  gap: 0.7rem;
  margin: 0.75rem 0.35rem 0;
  padding: 0.85rem;
  border: 1px solid rgb(20 125 131 / 16%);
  border-radius: 0.75rem;
  background:
    linear-gradient(135deg, rgb(20 125 131 / 7%), transparent),
    var(--bio-soft);
}

.drawer-status .status-dot {
  width: 0.5rem;
  height: 0.5rem;
  flex: 0 0 auto;
  border-radius: 50%;
  background: #49a36f;
  box-shadow: 0 0 0 0.25rem rgb(73 163 111 / 12%);
}

.drawer-status div {
  display: flex;
  min-width: 0;
  flex-direction: column;
}

.drawer-status strong {
  color: var(--bio-ink);
  font-size: 0.72rem;
}

.drawer-status small {
  margin-top: 0.2rem;
  color: var(--bio-muted);
  font-size: 0.65rem;
}

@media (max-width: 37.49rem) {
  .site-toolbar {
    min-height: 3.5rem;
    padding-inline: 0.75rem;
  }

  .brand {
    margin-left: 0.35rem;
  }

  .brand-copy small {
    display: none;
  }

  .brand-copy strong {
    font-size: 0.98rem;
  }
}
</style>