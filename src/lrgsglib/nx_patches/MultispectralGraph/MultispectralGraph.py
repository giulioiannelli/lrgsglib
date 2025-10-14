"""Implementation of the :class:`MultispectralGraph` SignedGraph subclass."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Mapping, Tuple

import networkx as nx

from ..common import *
from ..SignedGraph.SignedGraph import SignedGraph
from .generators import multiplicative_cascade_graph, multiplicative_cascade_probability_matrix


def _cascade_default_params() -> Dict[str, Any]:
    """Return defaults for the multiplicative cascade generator."""
    return {
        "p1": MSG_P1,
        "p2": MSG_P2,
        "p3": MSG_P3,
        "p4": MSG_P4,
        "fraction": MSG_FRACTION,
        "iterations": MSG_ITERATIONS,
        "probabilities": None,
    }


@dataclass(frozen=True)
class GeneratorSpec:
    """Container describing a multispectral generator."""

    graph_fn: Callable[..., nx.Graph]
    defaults_factory: Callable[[], Dict[str, Any]]
    path: str


GENERATOR_REGISTRY: Dict[str, GeneratorSpec] = {
    "multiplicative_cascade": GeneratorSpec(
        graph_fn=multiplicative_cascade_graph,
        defaults_factory=_cascade_default_params,
        path=MSG_PATHS.get("multiplicative_cascade", MSG_PATH),
    ),
}


class MultispectralGraph(SignedGraph):
    """A SignedGraph generated via multispectral graph generators."""

    def __init__(
        self,
        generator: Mapping[str, Any] | str | None = None,
        msg_type: str | None = None,
        msg_gen_params: Mapping[str, Any] | None = None,
        stdFnameSFFX: str = MSG_STDFN,
        sgpathn: str = MSG_SGPATH,
        only_const_mode: bool = MSG_ONLY_CONST_MODE,
        **kwargs,
    ) -> None:
        self._verify_pflip(kwargs.get("pflip", 0.0))

        legacy_param_keys = (
            "p1",
            "p2",
            "p3",
            "p4",
            "fraction",
            "iterations",
            "probabilities",
            "periodic",
        )
        legacy_overrides: Dict[str, Any] = {}
        for key in legacy_param_keys:
            if key in kwargs:
                legacy_overrides[key] = kwargs.pop(key)

        if "periodic" in legacy_overrides:
            raise ValueError("The 'periodic' flag is no longer supported for multispectral graphs.")

        if isinstance(generator, str):
            if msg_type is not None:
                raise ValueError("Specify the generator kind either via 'generator' string or 'msg_type', not both.")
            msg_type = generator
            generator = None

        generator_config = self._normalize_generator_config(
            generator=generator,
            msg_type=msg_type,
            msg_gen_params=msg_gen_params,
            overrides=legacy_overrides,
        )

        self.only_const_mode = only_const_mode
        self.msg_type = self._normalize_type(generator_config["kind"])
        self.generator_spec = self._resolve_generator_spec(self.msg_type)
        self.msg_gen_params = self._init_generator_params(generator_config["params"])
        self.__init_stdFname__(stdFnameSFFX)
        self.sgpathn = self._resolve_storage_path(sgpathn)

        self.probability_matrix = None
        self.seed_probabilities: Tuple[float, float, float, float] | Tuple[()] = tuple()
        self.iterations: int | None = None
        self.fraction: float | None = None
        self.syshape: Tuple[int, int] | None = None
        self.syshapePth = self.msg_type

        self._configure_metadata()

        if self.fraction is not None and not (0.0 < self.fraction <= 1.0):
            raise ValueError("fraction must lie within the interval (0, 1].")

        if not self.only_const_mode:
            self.__init_graph__()
        else:
            self.G = nx.Graph()

        super().__init__(self.G, **kwargs)

    def __init_stdFname__(self, suffix: str = "") -> None:
        base_name = f"{MSG_PHTABB}{self.msg_type}"
        self.std_fname = base_name + suffix if suffix else base_name

    def __init_graph__(self) -> None:
        params = dict(self.msg_gen_params)
        self.H = self.generator_spec.graph_fn(**params)
        self.G = nx.convert_node_labels_to_integers(self.H)

    def get_expected_num_nodes(self) -> int:
        if isinstance(self.syshape, tuple):
            return int(self.syshape[0] * self.syshape[1])
        raise ValueError("The multispectral graph shape is not defined.")

    def _normalize_type(self, msg_type: str) -> str:
        return msg_type.strip().lower()

    def _resolve_generator_spec(self, msg_type: str) -> GeneratorSpec:
        try:
            return GENERATOR_REGISTRY[msg_type]
        except KeyError as exc:
            available = ", ".join(sorted(GENERATOR_REGISTRY))
            raise ValueError(
                f"Unknown multispectral generator '{msg_type}'. "
                f"Available types: {available}."
            ) from exc

    def _normalize_generator_config(
        self,
        generator: Mapping[str, Any] | str | None,
        msg_type: str | None,
        msg_gen_params: Mapping[str, Any] | None,
        overrides: Mapping[str, Any],
    ) -> Dict[str, Any]:
        params: Dict[str, Any] = {}
        kind: str | None = None

        if generator is not None:
            if isinstance(generator, str):
                kind = generator
            else:
                if not isinstance(generator, Mapping):
                    raise TypeError("generator must be a mapping with a 'kind' entry.")
                generator_data = dict(generator)
                kind = generator_data.pop("kind", None) or generator_data.pop("type", None)
                params_from_config = generator_data.pop("params", None)
                if params_from_config is not None:
                    if not isinstance(params_from_config, Mapping):
                        raise TypeError("generator['params'] must be a mapping.")
                    params.update(params_from_config)
                for key, value in generator_data.items():
                    params[key] = value

        if kind is None:
            if msg_type is not None:
                kind = msg_type
            else:
                kind = MSG_DEFAULT_TYPE

        if not isinstance(kind, str):
            raise TypeError("Multispectral generator kind must be a string.")

        if msg_gen_params:
            params.update(msg_gen_params)

        if overrides:
            params.update(dict(overrides))

        return {"kind": kind, "params": params}

    def _init_generator_params(self, overrides: Mapping[str, Any]) -> Dict[str, Any]:
        params = self.generator_spec.defaults_factory()
        for key, value in overrides.items():
            params[key] = value
        return params

    def _resolve_storage_path(self, sgpathn: str | None) -> str:
        path_suffix = self.generator_spec.path or f"msg_{self.msg_type}"
        return pth_join(sgpathn, path_suffix) if sgpathn else path_suffix

    def _configure_metadata(self) -> None:
        if self.msg_type == "multiplicative_cascade":
            p1 = self.msg_gen_params.get("p1")
            p2 = self.msg_gen_params.get("p2")
            p3 = self.msg_gen_params.get("p3")
            p4 = self.msg_gen_params.get("p4")
            if None in (p1, p2, p3, p4):
                raise ValueError("Cascade generator requires probabilities p1, p2, p3, and p4.")

            self.seed_probabilities = (float(p1), float(p2), float(p3), float(p4))
            for idx, key in enumerate(("p1", "p2", "p3", "p4")):
                self.msg_gen_params[key] = self.seed_probabilities[idx]

            self.fraction = float(self.msg_gen_params.get("fraction", MSG_FRACTION))
            self.msg_gen_params["fraction"] = self.fraction

            self.iterations = int(self.msg_gen_params.get("iterations", MSG_ITERATIONS))
            self.msg_gen_params["iterations"] = self.iterations

            matrix = self.msg_gen_params.get("probabilities")
            if matrix is None:
                matrix = multiplicative_cascade_probability_matrix(
                    *self.seed_probabilities,
                    iterations=self.iterations,
                )
                self.msg_gen_params["probabilities"] = matrix

            self.probability_matrix = matrix
            if hasattr(matrix, "shape") and len(matrix.shape) == 2:
                self.syshape = (int(matrix.shape[0]), int(matrix.shape[1]))
                total_nodes = int(self.syshape[0] * self.syshape[1])
                self.syshapePth = f"{self.msg_type}_N={total_nodes}"
        else:
            self.syshapePth = self.msg_type


__all__ = ["MultispectralGraph"]
