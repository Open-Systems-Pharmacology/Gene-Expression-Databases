#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <tag> <asset-path> [release-name]"
  echo "Example: $0 v3.0.2 'PK-Sim DBs/Mouse/GENEDB_mouse_ADME_ONLY_BgeeRelease_15_2.expressionDB.tar.gz' 'OSP Expression DB v3.0.2'"
  exit 1
fi

TAG="$1"
ASSET_PATH="$2"
RELEASE_NAME="${3:-$TAG}"

if [[ ! -f "$ASSET_PATH" ]]; then
  echo "Error: asset file not found: $ASSET_PATH"
  exit 1
fi

if [[ -z "${GH_TOKEN:-}" ]]; then
  echo "Error: GH_TOKEN is not set. Export a GitHub token with repo scope first."
  exit 1
fi

if ! git config --get remote.origin.url >/dev/null 2>&1; then
  echo "Error: no git remote 'origin' configured."
  exit 1
fi

origin_url="$(git config --get remote.origin.url)"
repo=""

if [[ "$origin_url" =~ github.com[:/]([^/]+/[^/.]+)(\.git)?$ ]]; then
  repo="${BASH_REMATCH[1]}"
else
  echo "Error: could not parse GitHub repo from origin URL: $origin_url"
  exit 1
fi

api_base="https://api.github.com"
release_endpoint="$api_base/repos/$repo/releases"

auth_header="Authorization: Bearer $GH_TOKEN"
accept_header="Accept: application/vnd.github+json"

echo "Resolving release for tag: $TAG"
release_json="$(curl -sS -H "$auth_header" -H "$accept_header" "$release_endpoint/tags/$TAG")"

if echo "$release_json" | grep -q '"message": "Not Found"'; then
  echo "Release for tag $TAG not found. Creating draft release..."

  create_payload=$(cat <<EOF
{
  "tag_name": "$TAG",
  "name": "$RELEASE_NAME",
  "draft": true,
  "prerelease": false,
  "generate_release_notes": true
}
EOF
)

  release_json="$(curl -sS -X POST \
    -H "$auth_header" \
    -H "$accept_header" \
    -H "Content-Type: application/json" \
    "$release_endpoint" \
    -d "$create_payload")"
fi

release_id="$(echo "$release_json" | sed -n 's/.*"id": \([0-9][0-9]*\).*/\1/p' | head -1)"
upload_url_raw="$(echo "$release_json" | sed -n 's/.*"upload_url": "\([^"]*\)".*/\1/p' | head -1)"
upload_url="${upload_url_raw%\{*}"

if [[ -z "$release_id" || -z "$upload_url" ]]; then
  echo "Error: failed to parse release metadata from GitHub response."
  exit 1
fi

asset_name="$(basename "$ASSET_PATH")"

# Remove existing asset with same name, if present
assets_json="$(curl -sS -H "$auth_header" -H "$accept_header" "$release_endpoint/$release_id/assets")"
asset_id="$(echo "$assets_json" | tr '{' '\n' | grep -F "\"name\":\"$asset_name\"" | sed -n 's/.*"id":\([0-9][0-9]*\).*/\1/p' | head -1)"

if [[ -n "$asset_id" ]]; then
  echo "Deleting existing asset with same name: $asset_name (id=$asset_id)"
  curl -sS -X DELETE -H "$auth_header" -H "$accept_header" "$release_endpoint/assets/$asset_id" >/dev/null
fi

echo "Uploading asset: $asset_name"
curl -sS -X POST \
  -H "$auth_header" \
  -H "$accept_header" \
  -H "Content-Type: application/gzip" \
  --data-binary @"$ASSET_PATH" \
  "$upload_url?name=$asset_name" >/dev/null

echo "Done. Asset uploaded to release tag: $TAG"
echo "Release URL: https://github.com/$repo/releases/tag/$TAG"
