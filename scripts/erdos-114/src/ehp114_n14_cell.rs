use serde::Serialize;
use serde_json::Value;

pub const DEFAULT_SUB_I: usize = 6;
pub const DEFAULT_SUB_J: usize = 4;
pub const ROOT_AFFINE_SUBDIVISION: usize = 8;
pub const U0_PARENT: [f64; 2] = [-0.0017499999999999998, 0.0];
pub const U1_PARENT: [f64; 2] = [0.0, 0.0017499999999999998];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct CellSpec {
    pub sub_i: usize,
    pub sub_j: usize,
}

impl CellSpec {
    pub fn new(sub_i: usize, sub_j: usize) -> Result<Self, String> {
        if sub_i >= ROOT_AFFINE_SUBDIVISION || sub_j >= ROOT_AFFINE_SUBDIVISION {
            return Err(format!(
                "invalid subcell ({sub_i},{sub_j}); expected 0 <= sub_i, sub_j < {ROOT_AFFINE_SUBDIVISION}"
            ));
        }
        Ok(Self { sub_i, sub_j })
    }

    pub fn hard_cell() -> Self {
        Self {
            sub_i: DEFAULT_SUB_I,
            sub_j: DEFAULT_SUB_J,
        }
    }

    pub fn is_hard_cell(self) -> bool {
        self.sub_i == DEFAULT_SUB_I && self.sub_j == DEFAULT_SUB_J
    }

    pub fn tag(self) -> String {
        format!("CELL-{:02}-{:02}", self.sub_i, self.sub_j)
    }

    pub fn u0_interval(self) -> [f64; 2] {
        split_interval_index(U0_PARENT, ROOT_AFFINE_SUBDIVISION, self.sub_i)
    }

    pub fn u1_interval(self) -> [f64; 2] {
        split_interval_index(U1_PARENT, ROOT_AFFINE_SUBDIVISION, self.sub_j)
    }
}

pub fn split_interval_index(pair: [f64; 2], n: usize, index: usize) -> [f64; 2] {
    let step = (pair[1] - pair[0]) / (n as f64);
    [
        pair[0] + (index as f64) * step,
        pair[0] + ((index + 1) as f64) * step,
    ]
}

fn parse_usize_value(value: &Value, label: &str) -> Result<usize, String> {
    value
        .as_u64()
        .map(|n| n as usize)
        .ok_or_else(|| format!("subcell field {label} is not a nonnegative integer"))
}

pub fn source_subcell_from_value(value: &Value) -> Result<Option<CellSpec>, String> {
    let Some(subcell) = value.get("subcell") else {
        return Ok(None);
    };
    if let Some(arr) = subcell.as_array() {
        if arr.len() != 2 {
            return Err("source subcell array does not have length 2".to_string());
        }
        return CellSpec::new(
            parse_usize_value(&arr[0], "subcell[0]")?,
            parse_usize_value(&arr[1], "subcell[1]")?,
        )
        .map(Some);
    }
    if let Some(obj) = subcell.as_object() {
        let sub_i = obj
            .get("sub_i")
            .ok_or_else(|| "source subcell object missing sub_i".to_string())
            .and_then(|v| parse_usize_value(v, "sub_i"))?;
        let sub_j = obj
            .get("sub_j")
            .ok_or_else(|| "source subcell object missing sub_j".to_string())
            .and_then(|v| parse_usize_value(v, "sub_j"))?;
        return CellSpec::new(sub_i, sub_j).map(Some);
    }
    Err("source subcell is neither [i,j] nor {sub_i,sub_j}".to_string())
}

pub fn ensure_source_subcell_matches(
    source: &Value,
    requested: CellSpec,
) -> Result<String, String> {
    match source_subcell_from_value(source)? {
        Some(source_cell) if source_cell == requested => Ok("SOURCE_SUBCELL_MATCH".to_string()),
        Some(source_cell) => Err(format!(
            "SOURCE_SUBCELL_MISMATCH: requested {} but source artifact declares {}",
            requested.tag(),
            source_cell.tag()
        )),
        None if requested.is_hard_cell() => {
            Ok("SOURCE_SUBCELL_LEGACY_HARD_CELL_ASSUMED".to_string())
        }
        None => Err(format!(
            "SOURCE_SUBCELL_UNVERIFIED: requested {} but source artifact lacks subcell metadata",
            requested.tag()
        )),
    }
}
