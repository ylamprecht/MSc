-- ============================================================
-- Export the curated literature variant set
--
-- For each variant, one odds ratio is selected from the papers
-- reporting it, ranked by evidence level, then sample size,
-- then p-value.
--
-- Output: published_risk_variants.csv
-- ============================================================

WITH ranked_or AS (
    SELECT
        sv.RSID,
        vp.Odds_ratio,
        vp.p_value,
        p.Paper_ID,
        p.DOI,
        p.Evidence_level,
        pop.Sample_size,

        ROW_NUMBER() OVER (
            PARTITION BY sv.RSID
            ORDER BY
                -- 1. Highest evidence level
                CASE p.Evidence_level
                    WHEN 'High' THEN 1
                    WHEN 'Moderate' THEN 2
                    WHEN 'Candidate' THEN 3
                    ELSE 4
                END,

                -- 2. Prefer papers with a reported sample size
                CASE
                    WHEN pop.Sample_size IS NULL THEN 1
                    ELSE 0
                END,

                -- 3. Within the same evidence level,
                -- choose the largest sample size
                pop.Sample_size DESC,

                -- 4. If sample size is unavailable,
                -- choose the smallest p-value
                CASE
                    WHEN pop.Sample_size IS NULL
                    THEN vp.p_value
                    ELSE NULL
                END ASC,

                -- Final tie-breaker
                p.Paper_ID ASC
        ) AS rn

    FROM Significant_variants AS sv

    JOIN Variants_Papers AS vp
        ON sv.Significant_variant_ID = vp.Significant_variant_ID

    JOIN Paper AS p
        ON vp.Paper_ID = p.Paper_ID

    LEFT JOIN Population AS pop
        ON p.Paper_ID = pop.Paper_ID

    -- Only rank papers that report an OR
    WHERE vp.Odds_ratio IS NOT NULL
)

SELECT
    vd.RSID,
    vd.Gene,
    vd.Chromosome,
    vd.Position,
    vd.Mutation_type,
    vd.Reference_allele,
    vd.Associated_allele,
    sv.Association,
    sv.Disease_status,
    ro.Odds_ratio,
    ro.p_value,
    ro.DOI

FROM Variant_detail AS vd

JOIN Significant_variants AS sv
    ON vd.RSID = sv.RSID

LEFT JOIN ranked_or AS ro
    ON vd.RSID = ro.RSID
    AND ro.rn = 1

ORDER BY
    vd.Chromosome,
    vd.Position;
