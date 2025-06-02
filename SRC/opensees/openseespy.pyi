from typing import overload, Literal

class _Elements:

    @overload
    def element(self, type: str, tag: int, nodes: tuple, **kwargs) -> int:
        """
        Create an element of the specified type.
        """

    @overload
    def element(self,
                __elmType: Literal["PrismFrame"],
                tag: int,
                nodes: tuple,
                section: int,
                transform: int,
                *args
            ) -> int:
        """
        Create a linear elastic prismatic frame element.        
        
        :param tag: unique :ref:`element` tag
        :type tag: |integer|
        :param nodes: tuple of *two* integer :ref:`node` tags
        :type nodes: tuple
        :param section: tag of a previously-defined :ref:`section`
        :type section: |integer|
        :param transform: identifier for previously-defined coordinate-transformation
        :type transform: |integer|
        """

    @overload
    def element(self,
                __elmType: Literal["forceBeamColumn"],
                tag: int,
                nodes: tuple,
                nip: int,
                transform: int,
                *args
               ) -> int: ...



class Model(_Elements):
    pass
