# Attempt the two declared BL21 guides independently.
# A non-zero runner status records a failure and still invokes the other guide.
attempt_peer_guides() {
    local runner="$1"
    local guide_id
    ATTEMPTED_GUIDES=""
    GUIDE_FAILURES=0
    for guide_id in BL21_G1 BL21_G2; do
        if [[ -n "$ATTEMPTED_GUIDES" ]]; then
            ATTEMPTED_GUIDES+=" "
        fi
        ATTEMPTED_GUIDES+="$guide_id"
        if ! "$runner" "$guide_id"; then
            GUIDE_FAILURES=$((GUIDE_FAILURES + 1))
            continue
        fi
    done
}
