#!/usr/bin/env bash

set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly GAPS_REPOSITORY="https://github.com/hongliangduan/GAPS.git"
readonly GAPS_REVISION="c67cba6ca3ef45470bee8bbb2f0c832680c18044"
readonly PEPGLAD_REPOSITORY="https://github.com/THUNLP-MT/PepGLAD.git"
readonly PEPGLAD_REVISION="bad015ca50c312a89482adb5220c3d907f13df5c"
readonly CHECKPOINT_URL="https://github.com/THUNLP-MT/PepGLAD/releases/download/v1.0/checkpoints.zip"
readonly CHECKPOINT_SHA256="a610e079492a1b1ceab213aea3a0ab875846415e16470897bb54850557f462a0"

checkpoint_archive=""

cleanup() {
    if [[ -n "${checkpoint_archive}" && -f "${checkpoint_archive}" ]]; then
        rm -f -- "${checkpoint_archive}"
    fi
}
trap cleanup EXIT

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

install_repository() {
    local name="$1"
    local repository="$2"
    local revision="$3"
    local destination="${SCRIPT_DIR}/${name}"
    local installed_revision

    if [[ ! -e "${destination}" ]]; then
        printf 'Installing %s...\n' "${name}"
        git clone --quiet "${repository}" "${destination}"
        git -C "${destination}" checkout --quiet --detach "${revision}"
        return
    fi

    [[ -d "${destination}/.git" ]] ||
        die "${destination} exists but is not a Git repository"

    installed_revision="$(git -C "${destination}" rev-parse HEAD)"
    [[ "${installed_revision}" == "${revision}" ]] ||
        die "${name} is at ${installed_revision}; expected ${revision}. Move or remove ${destination}, then rerun this installer."

    printf 'Using existing %s checkout.\n' "${name}"
}

apply_gaps_patch() {
    local gaps_directory="${SCRIPT_DIR}/GAPS"
    local patch_file="${SCRIPT_DIR}/patches/GAPS.patch"

    [[ -f "${patch_file}" ]] || die "GAPS patch not found: ${patch_file}"
    touch "${gaps_directory}/__init__.py"

    if git -C "${gaps_directory}" apply --reverse --check \
        --ignore-space-change "${patch_file}" >/dev/null 2>&1; then
        printf 'GAPS compatibility patch is already applied.\n'
    elif git -C "${gaps_directory}" apply --check \
        --ignore-space-change "${patch_file}" >/dev/null 2>&1; then
        git -C "${gaps_directory}" apply --ignore-space-change "${patch_file}"
        printf 'Applied the GAPS compatibility patch.\n'
    else
        die "The GAPS compatibility patch cannot be applied cleanly"
    fi
}

download_file() {
    local url="$1"
    local destination="$2"

    if command -v curl >/dev/null 2>&1; then
        curl --fail --location --retry 3 --show-error \
            --output "${destination}" "${url}"
    else
        wget --tries=3 --output-document="${destination}" "${url}"
    fi
}

sha256_digest() {
    local file="$1"
    local digest
    local ignored

    if command -v sha256sum >/dev/null 2>&1; then
        read -r digest ignored < <(sha256sum "${file}")
    else
        read -r digest ignored < <(shasum -a 256 "${file}")
    fi
    printf '%s' "${digest}"
}

install_pepglad_checkpoints() {
    local checkpoints_directory="${SCRIPT_DIR}/PepGLAD/checkpoints"
    local checkpoint
    local actual_checksum
    local missing_checkpoint=false
    local required_checkpoints=(
        codesign.ckpt
        codesign_pepbdb.ckpt
        fixseq.ckpt
        fixseq_pepbdb.ckpt
    )

    for checkpoint in "${required_checkpoints[@]}"; do
        if [[ ! -s "${checkpoints_directory}/${checkpoint}" ]]; then
            missing_checkpoint=true
            break
        fi
    done

    if [[ "${missing_checkpoint}" == false ]]; then
        printf 'PepGLAD checkpoints are already installed.\n'
        return
    fi

    printf 'Downloading PepGLAD v1.0 checkpoints...\n'
    checkpoint_archive="$(mktemp /tmp/allopep-checkpoints.XXXXXX.zip)"
    download_file "${CHECKPOINT_URL}" "${checkpoint_archive}"

    actual_checksum="$(sha256_digest "${checkpoint_archive}")"
    [[ "${actual_checksum}" == "${CHECKPOINT_SHA256}" ]] ||
        die "PepGLAD checkpoint checksum mismatch (got ${actual_checksum})"

    unzip -oq "${checkpoint_archive}" 'checkpoints/*' \
        -d "${SCRIPT_DIR}/PepGLAD"

    for checkpoint in "${required_checkpoints[@]}"; do
        [[ -s "${checkpoints_directory}/${checkpoint}" ]] ||
            die "Checkpoint was not extracted: ${checkpoint}"
    done

    rm -f -- "${checkpoint_archive}"
    checkpoint_archive=""
    printf 'Installed PepGLAD checkpoints.\n'
}

main() {
    require_command git
    require_command unzip
    if ! command -v curl >/dev/null 2>&1 &&
        ! command -v wget >/dev/null 2>&1; then
        die 'Required command not found: install curl or wget'
    fi
    if ! command -v sha256sum >/dev/null 2>&1 &&
        ! command -v shasum >/dev/null 2>&1; then
        die 'Required command not found: install sha256sum or shasum'
    fi

    install_repository GAPS "${GAPS_REPOSITORY}" "${GAPS_REVISION}"
    apply_gaps_patch
    install_repository PepGLAD "${PEPGLAD_REPOSITORY}" "${PEPGLAD_REVISION}"
    install_pepglad_checkpoints

    printf 'AlloPep dependencies installed successfully.\n'
}

main "$@"
